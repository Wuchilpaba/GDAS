import pandas as pd
import numpy as np
import ray
import pickle
import pymc as pm
import arviz as az
import pathlib, os

# Threshold for deciding methods
BOOTSTRAP_THRESHOLD = 20  # minimum samples per group for bootstrap

def compute_cohens_d(means, variances, n1, n2):
    """
    Calculate Cohen's d effect size between two groups.
    """
    observed_diff = means.iloc[1] - means.iloc[0]
    pooled_var = ((n1 - 1) * variances.iloc[0] + (n2 - 1) * variances.iloc[1]) / max(n1 + n2 - 2, 1)
    pooled_std = np.sqrt(pooled_var)
    return observed_diff / pooled_std if pooled_std != 0 else np.nan


def compute_hedges_g(cohen_d, n1, n2):
    """
    Apply Hedges' g correction to Cohen's d.
    """
    denom = 4 * (n1 + n2) - 9
    if denom <= 0:
        return np.nan
    correction = 1 - 3 / denom
    return cohen_d * correction


def compute_effect_size(means, variances, sizes):
    """
    Always compute Cohen's d, then optionally Hedges' g correction.
    The decision to use bootstrap vs Bayesian is made separately based on total samples.
    """
    groups = sorted(means.index.tolist())
    n1, n2 = sizes[groups[0]], sizes[groups[1]]
    d = compute_cohens_d(means.loc[groups], variances.loc[groups], n1, n2)
    return d


def robust_bayesian_analysis(
        data: pd.DataFrame,
        alpha: float = 0.05,
        tune: int = 2000,
        draws: int = 2000,
        min_effective_samples: int = 500
) -> tuple[float, float, float, float]:
    """
    鲁棒的贝叶斯差异分析函数，适用于极端小样本场景

    参数：
    data : pd.DataFrame
        必须包含'Group'和'Intensity'两列
    alpha : float, 可选
        显著性水平 (默认0.05)
    tune : int, 可选
        MCMC预热迭代次数 (默认2000)
    draws : int, 可选
        正式采样次数 (默认2000)
    min_effective_samples : int, 可选
        最小有效样本量阈值 (默认500)

    返回：
    Tuple[float, float, float, float]
        (效应量均值, HDI下限, HDI上限, 贝叶斯p值)
        无效情况返回四个np.nan
    """
    # ===================== 输入验证 =====================
    try:
        if not {'Group', 'Intensity'}.issubset(data.columns):
            raise ValueError("输入数据必须包含'Group'和'Intensity'列")

        groups = data['Group'].unique()
        if len(groups) != 2:
            raise ValueError("需要且仅需要两个实验组")

        # ================= 数据预处理 ==================
        # 提取数据并处理缺失值
        obs0_raw = data.loc[data['Group'] == groups[0], 'Intensity'].dropna().values
        obs1_raw = data.loc[data['Group'] == groups[1], 'Intensity'].dropna().values

        # Tukey异常值过滤（内置实现）
        obs0 = _tukey_filter(obs0_raw)
        obs1 = _tukey_filter(obs1_raw)

        # ============== 自适应先验设置 ================
        combined = np.concatenate([obs0, obs1])
        n_obs = len(combined)

        # 使用稳健的统计量计算
        mu_center = np.mean(combined) if n_obs > 0 else 0.0
        mu_sd = _calc_mu_sd(combined, n_obs)
        sigma_prior = _calc_sigma_prior(combined, n_obs)

        # =============== 模型构建 ================
        with pm.Model() as robust_model:
            # 使用Cauchy分布增强鲁棒性
            mu0 = pm.Cauchy('mu0', alpha=mu_center, beta=mu_sd)
            mu1 = pm.Cauchy('mu1', alpha=mu_center, beta=mu_sd)
            sigma = pm.HalfCauchy('sigma', beta=sigma_prior)

            # 观测值使用T分布
            if len(obs0) > 0:
                pm.StudentT('obs0', nu=4, mu=mu0, sigma=sigma, observed=obs0)
            if len(obs1) > 0:
                pm.StudentT('obs1', nu=4, mu=mu1, sigma=sigma, observed=obs1)

            # 效应量计算增加稳定性保护
            pooled_sigma = pm.math.sqrt((sigma**2 + sigma**2) / 2 + 1e-9)
            cohen_d = pm.Deterministic('cohen_d', (mu1 - mu0) / pooled_sigma)

            # ============ 自适应采样策略 ============
            trace = pm.sample(
                tune=3000,
                draws=3000,
                chains=4,
                cores=2,
                target_accept=0.95,
                return_inferencedata=True
            )

            # ============ 收敛诊断 ============
            if 'cohen_d' not in trace.posterior:
                raise RuntimeError("效应量变量未找到")

            ess = az.ess(trace, var_names=['cohen_d'])['cohen_d'].values.min()
            print(az.ess)
            print(type(az.ess))
            rhat = az.rhat(trace, var_names=['cohen_d'])['cohen_d'].values.max()
            print(az.rhat)
            print(type(az.rhat))
            if ess < min_effective_samples:
                raise RuntimeError(f"有效样本量不足: {ess} < {min_effective_samples}")
            if rhat > 1.05:
                raise RuntimeError(f"未收敛, R-hat值过高: {rhat:.3f}")

        # ============ 结果计算 ============
        d_samples = trace.posterior['cohen_d'].values.flatten()
        d_samples = d_samples[~np.isnan(d_samples)]

        if len(d_samples) < min_effective_samples:
            raise RuntimeError("有效样本不足")

        hdi = az.hdi(d_samples, hdi_prob=1 - alpha)
        p_value = 2 * min((d_samples < 0).mean(), (d_samples > 0).mean())

        return (
            float(np.mean(d_samples)),
            float(hdi[0]) if not np.isnan(hdi).all() else np.nan,
            float(hdi[1]) if not np.isnan(hdi).all() else np.nan,
            float(p_value)
        )

    except KeyError as ke:
        print(f"数据列错误: {str(ke)}")
    except ValueError as ve:
        print(f"输入验证失败: {str(ve)}")
    except RuntimeError as re:
        print(f"运行时错误: {str(re)}")
    except Exception as e:
        print(f"未捕获异常: {str(e)}")

    return (np.nan, np.nan, np.nan, np.nan)


def _tukey_filter(data: np.ndarray, k: float = 1.5) -> np.ndarray:
    """Tukey异常值过滤"""
    if len(data) < 4:
        return data  # 样本过少不处理

    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower = q1 - k * iqr
    upper = q3 + k * iqr
    return data[(data >= lower) & (data <= upper)]


def _calc_mu_sd(combined: np.ndarray, n_obs: int) -> float:
    """动态计算均值先验标准差"""
    if n_obs == 0:
        return 5.0  # 默认先验
    elif n_obs == 1:
        return max(np.ptp(combined), 1.0) if len(combined) > 0 else 5.0
    else:
        return 10 * np.std(combined)


def _calc_sigma_prior(combined: np.ndarray, n_obs: int) -> float:
    """动态计算方差先验"""
    if n_obs <= 1:
        return 5.0
    else:
        return max(np.std(combined), 0.5)

@ray.remote
def bootstrap_analysis(entries: pd.DataFrame) -> tuple:
    """
    Perform bootstrap-based effect size calculation and its permutation.
    Returns a tuple with ('bootstrap', eff, nan, nan, nan, perm_eff).
    """
    df = entries.copy()
    # 统计每个组的大小和组别
    sizes = {g: len(df[df['Group'] == g]) for g in df['Group'].unique()}
    groups = sorted(sizes.keys())

    # 生成 bootstrap 样本并计算均值和方差
    boot = pd.concat([
        df[df['Group'] == g].sample(n=sizes[g], replace=True)
        for g in groups
    ])
    bm = boot.groupby('Group')['Intensity'].mean()
    bv = boot.groupby('Group')['Intensity'].var()
    eff = compute_effect_size(bm, bv, sizes)

    # 置换检验计算 effect size
    perm_labels = np.random.permutation(df['Group'].values)
    perm_df = df.copy()
    perm_df['Group'] = perm_labels
    pmn = perm_df.groupby('Group')['Intensity'].mean()
    pvv = perm_df.groupby('Group')['Intensity'].var()
    ps = {g: int((perm_labels == g).sum()) for g in groups}
    perm_eff = compute_effect_size(pmn, pvv, ps)

    return ('bootstrap', eff, np.nan, np.nan, np.nan, perm_eff)


def analyze_iteration(entries: pd.DataFrame, n_iter: int = 1000):
    """
    根据总样本量决定使用 Bootstrap 还是贝叶斯分析。
    若总样本数 >= BOOTSTRAP_THRESHOLD，则调用远程的 bootstrap_analysis；
    否则执行本地鲁棒贝叶斯分析。
    """
    df = entries.copy()
    # 检查组别数量
    sizes = {g: len(df[df['Group'] == g]) for g in df['Group'].unique()}
    if len(sizes) < 2:
        return ('invalid', np.nan, np.nan, np.nan, np.nan, np.nan)

    total_samples = len(df)
    # 根据总样本数选择方法
    if total_samples >= BOOTSTRAP_THRESHOLD:
        # 使用 Ray 远程执行 bootstrap
        future = [bootstrap_analysis.remote(entries) for _ in range(n_iter)]
        results = ray.get(future)
        return results
    else:
        # 使用鲁棒贝叶斯分析
        eff, low, high, bp = robust_bayesian_analysis(df)
        return [('bayesian', eff, low, high, bp, np.nan)]

def main():
    ray.init(ignore_reinit_error=True)
    try:
        with open('data_cache_2.pickle','rb') as f:
            data = pickle.load(f)
    except FileNotFoundError:
        print("Error: Data file not found.")
        return

    results = analyze_iteration(data,n_iter=1000)

    out_path = pathlib.Path(os.getcwd()) / "result_cache_2.pickle"
    with open(out_path, 'wb') as f:
        pickle.dump(results, f)

    print(f"Total iterations: {len(results)}")
    # bcount = sum(1 for r in results if r[0]=='bootstrap')
    # baycount = sum(1 for r in results if r[0]=='bayesian')
    # print(f"Bootstrap iterations: {bcount}")
    # print(f"Bayesian iterations: {baycount}")
    sample = next(r for r in results if r[0]=='bayesian')
    print(f"Example Bayesian HDI: [{sample[2]}, {sample[3]}]")
    ray.shutdown()

if __name__=='__main__':
    main()
