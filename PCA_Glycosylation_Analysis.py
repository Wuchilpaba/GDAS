import pickle, pathlib
import time

from Bio import ExPASy, SwissProt
from scipy import stats
from xgboost import XGBRegressor

from resource_path import resource_path
import ctypes
# import modin.config as cfg
import pandas as pd
import numpy as np
import re, os
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from AnalysisSignalModule import AnalysisSignal
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial import distance
from scipy.stats import chi2
from sklearn.metrics import mean_squared_error, r2_score, silhouette_score
from sklearn.model_selection import cross_val_score, GridSearchCV, StratifiedKFold
import Annotation
import docker_ray_manager_1, docker_ray_manager_2

path = pathlib.PurePosixPath('/usr/bin/python')
final_signal = AnalysisSignal()
mpl.rcParams['font.family'] = 'Times New Roman'


def compute_cohens_d(means, variances, group_sizes):
    groups = sorted(means.index.tolist())
    n1 = group_sizes.iloc[0] if isinstance(group_sizes, pd.Series) else group_sizes.get(groups[0], 0)
    n2 = group_sizes.iloc[1] if isinstance(group_sizes, pd.Series) else group_sizes.get(groups[1], 0)
    var1 = variances.iloc[0] if isinstance(variances, pd.Series) else variances.get(groups[0], 0)
    var2 = variances.iloc[1] if isinstance(variances, pd.Series) else variances.get(groups[1], 0)
    # 检查分母有效性
    if (n1 + n2) <= 2:
        return np.nan
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    d = (means[groups[1]] - means[groups[0]]) / np.sqrt(pooled_var) if pooled_var > 0 else np.nan
    return d


def compute_hedges_g(cohen_d, n1, n2):
    denominator = 4 * (n1 + n2) - 9
    if denominator <= 0:
        return np.nan
    correction = 1 - 3 / denominator
    return cohen_d * correction


def compute_effect_size(means, variances, group_sizes):
    groups = sorted(means.index.tolist())
    if len(groups) != 2:
        return np.nan, "Invalid groups"
    n1, n2 = group_sizes.get(groups[0], 0), group_sizes.get(groups[1], 0)
    cohen_d = compute_cohens_d(means, variances, {groups[0]: n1, groups[1]: n2})
    if min(n1, n2) < 20 and not np.isnan(cohen_d):
        return compute_hedges_g(cohen_d, n1, n2), "Hedge's g"
    return cohen_d, "Cohen's d"


def calculate_pvalue_and_ci(observed, bootstrap_diffs, original_data, alpha=0.05):
    """计算双尾p值及BCa校正置信区间"""
    p_value = np.mean(np.abs(bootstrap_diffs) >= np.abs(observed))

    n_samples = len(original_data)

    # 计算偏差校正因子z0（避免极端值）
    mean_less = np.mean(bootstrap_diffs < observed)
    epsilon = 1e-10  # 微小偏移量防止z0为无穷大
    mean_less = np.clip(mean_less, epsilon, 1 - epsilon)
    z0 = stats.norm.ppf(mean_less)

    # 计算Jackknife效应量
    jack_effects = []
    for i in range(n_samples):
        row_label = original_data.index[i]
        jack_sample = original_data.drop(index=row_label)
        jack_means = jack_sample.groupby('Group')['Intensity'].mean()
        jack_vars = jack_sample.groupby('Group')['Intensity'].var()
        jack_sizes = jack_sample.groupby('Group')['Intensity'].count()
        effect, _ = compute_effect_size(jack_means, jack_vars, jack_sizes)
        jack_effects.append(effect)

    # 计算加速因子a（避免除以零）
    mean_jack = np.mean(jack_effects)
    jn = np.array(jack_effects) - mean_jack
    sum_jn2 = np.sum(jn ** 2)
    if sum_jn2 == 0:
        # 所有效应量相同，加速因子a=0，BCa退化为普通百分位数
        a = 0.0
    else:
        a = np.sum(jn ** 3) / (6 * sum_jn2 ** 1.5)

    # 计算调整后的分位数（带保护条件）
    z_alpha = stats.norm.ppf([alpha / 2, 1 - alpha / 2])
    denominator = 1 - a * (z0 + z_alpha)

    if np.any(denominator <= 0) or np.any(np.isnan(denominator)):
        # 分母无效，回退到普通百分位数法
        lower = np.percentile(bootstrap_diffs, alpha / 2 * 100)
        upper = np.percentile(bootstrap_diffs, (1 - alpha / 2) * 100)
    else:
        alpha_vals = stats.norm.cdf(z0 + (z0 + z_alpha) / denominator)
        alpha_vals = np.clip(alpha_vals, 0.0, 1.0)
        # 再次检查alpha_vals有效性
        if np.any(np.isnan(alpha_vals)):
            lower = np.percentile(bootstrap_diffs, alpha / 2 * 100)
            upper = np.percentile(bootstrap_diffs, (1 - alpha / 2) * 100)
        else:
            lower = np.percentile(bootstrap_diffs, alpha_vals[0] * 100)
            upper = np.percentile(bootstrap_diffs, alpha_vals[1] * 100)

    return p_value, (lower, upper)


def connection_docker_1(data):
    worker_ray_path = pathlib.Path(resource_path('resources'))
    with open(f'{worker_ray_path}\\data_cache1.pickle', 'wb') as f:
        pickle.dump(data, f)


def connection_docker_2(data):
    worker_ray_path = pathlib.Path(resource_path('resources'))
    with open(f'{worker_ray_path}\\data_cache_2.pickle', 'wb') as f:
        pickle.dump(data, f)


class byonic_r:
    def __init__(self, file_list1, file_list2):
        self.file_list1 = file_list1
        self.file_list2 = file_list2
        self.data = self.byonic_multiple_data()

    def read_byonic_selection(self, p):
        pattern = r'\|([A-Z0-9]+)\|'
        # ray.init(local_mode=True)
        # 读取Excel文件中的Spectra和Proteins工作表
        df = pd.read_excel(p, engine='openpyxl', sheet_name='Spectra')
        df_P = pd.read_excel(p, engine='openpyxl', sheet_name='Proteins')
        final_signal.write(p)
        # df = df_all['Spectra']
        # df_P = df_all['Protein']
        # 提取Proteins工作表中的Description和Total Intensity
        P = list(filter(lambda pro: '>Reverse' not in pro, list(df_P['Description'])))
        # filter(lambda pro: '>Reverse' not in pro, list(df_P)
        # 使用正则表达式进行替换
        updated_Protein = []
        for Pro in P:
            match = re.search(pattern, Pro)
            if match:
                # 如果匹配到，提取中间部分如 O15031
                updated_Protein.append(match.group(1))
            else:
                # 如果没有匹配到，保留原字符串
                updated_Protein.append(Pro)
        P = updated_Protein
        I = list(df_P["Total\nIntensity"])
        # 使用正则表达式筛选Spectra中的Glycan列
        pattern = r'^Glycan.*'
        columns = df.columns
        filtered_columns = [col for col in columns if re.match(pattern, col)]

        # 根据筛选后的列名创建新的 DataFrame
        filtered_df = df[df[filtered_columns[0]].notnull()]
        filtered_df = filtered_df[~filtered_df['Protein Name'].str.contains('>Reverse')]
        # 提取Spectra中的Protein Name和Glycan列
        Glycan = list(filtered_df[filtered_columns[0]])
        Protein = list(filtered_df['Protein Name'])
        # 使用正则表达式进行替换
        updated_Protein = []
        pattern = r'\|([A-Z0-9]+)\|'
        for Pro in Protein:
            match = re.search(pattern, Pro)
            if match:
                # 如果匹配到，提取中间部分如 O15031
                updated_Protein.append(match.group(1))
            else:
                # 如果没有匹配到，保留原字符串
                updated_Protein.append(Pro)
        Protein = updated_Protein
        Proteins = set(Protein)

        # 创建反射字典：将蛋白质名称映射到强度
        reflexion_P = {}
        for index, Pr in enumerate(P):
            reflexion_P[Pr] = I[index]
            if Pr not in Proteins:
                del reflexion_P[Pr]

        # 创建字典：将蛋白质名称映射到Glycan
        reflexion_G = {}
        element_count = Counter(Protein)
        for i, element in enumerate(Protein):
            if element_count[element] > 1:
                if element not in reflexion_G:
                    reflexion_G[element] = []
                reflexion_G[element].append(Glycan[i])
            else:
                reflexion_G[element] = Glycan[i]
        self.protein_list = reflexion_P.keys()

        # 返回糖蛋白与强度（reflexion_P）和糖蛋白与聚糖（reflexion_G）的字典
        return reflexion_P, reflexion_G, reflexion_P.keys()

    def byonic_multiple_data(self):
        group_li = [self.file_list1, self.file_list2]
        combined_data_final = []

        # 存储所有文件数据
        for index, item in enumerate(group_li):
            data = []
            # 遍历每个组别的文件
            for file in item:
                path = pathlib.Path(file)
                reflexion_P, reflexion_G, protein_list = self.read_byonic_selection(path)  # 调用read_byonic_selection函数

                # 创建新的DataFrame，将Proteins和Intensity结合，并添加组别信息
                df = pd.DataFrame({
                    'Protein': list(reflexion_P.keys()),
                    'Intensity': list(reflexion_P.values()),
                    'Glycan': [reflexion_G.get(p, None) for p in reflexion_P.keys()],
                    'Group': index + 1  # 1表示实验组，2表示对照组
                })
                data.append(df)

            # 合并每个组别的所有文件数据
            combined_data = pd.concat(data, ignore_index=True)
            combined_data_final.append(combined_data)

        # 合并所有组别的数据为一个完整的DataFrame
        final_data = pd.concat(combined_data_final, ignore_index=True)
        return final_data

    def bootstrap_analysis(self):
        # 存储每个糖蛋白的Bootstrap结果
        all_results = {}

        # 按Protein分组
        for protein, group_data in self.data.groupby('Protein'):
            results = {}

            # 获取组别
            groups = group_data['Group'].unique()
            if len(groups) != 2:
                # raise ValueError("Bootstrap analysis requires exactly two groups.")
                all_results[protein] = None
            else:
                # 分离两组数据
                group1_data = group_data[group_data['Group'] == groups[0]]['Intensity'].values
                group2_data = group_data[group_data['Group'] == groups[1]]['Intensity'].values

                # 使用StandardScaler标准化每组数据
                scaler = StandardScaler()
                # Reshape for StandardScaler (it expects 2D array)
                group1_data_scaled = scaler.fit_transform(group1_data.reshape(-1, 1)).flatten()
                group2_data_scaled = scaler.fit_transform(group2_data.reshape(-1, 1)).flatten()

                data = (group1_data_scaled, group2_data_scaled)

                # Launch n_iterations bootstrap iterations in parallel

                final_signal.write('Deploying bootstrap iteration tasks on Ray cluster....')
                connection_docker_1(data)
                container = docker_ray_manager_1.start_ray_container()
                current_directory = os.getcwd()
                worker_ray_path = pathlib.Path(resource_path("resources")) / "result_cache1.pickle"
                # log = pathlib.Path(current_directory,"docker_ray",'task2_complete.txt')
                while not os.path.exists(worker_ray_path):
                    time.sleep(3)
                final_signal.write(f'{protein} task1 has finished')
                with open(worker_ray_path, 'rb') as handle:
                    boot_means_1 = pickle.load(handle)
                # Convert the list of results to a numpy array
                boot_means = np.array(boot_means_1)
                docker_ray_manager_1.stop_ray_container(container)

                # 计算置信区间
                conf_interval = np.percentile(boot_means, [2.5, 97.5])
                mean_diff = np.mean(boot_means)

                # 计算P值
                # p_value = np.mean(np.abs(boot_means) >= np.abs(original_mean_diff))

                # 存储结果
                results['Mean Difference'] = mean_diff
                results['95% CI'] = conf_interval
                # results['P-Value'] = p_value  # 添加P值

                # 更新结果字典
                all_results[protein] = results

        # 返回每个糖蛋白的Bootstrap结果
        return all_results


# 输入一种糖蛋白的Glycan-Intensity Ratio，Sites-Intensity Ratio
# 或许字典更为合适
# Protein_list comes from Byonic Glycosylation list
# glycan_list comes from Analysis module of GlycReSoft/O-Pair Result
class Vector:
    def __init__(self, glyc: str, expression_level_path: str, protein_list, glycan_list):
        self.glycan = glyc[0].replace('/', '\\')
        self.expression_level = expression_level_path[0].replace('/', '\\')
        self.protein_list = protein_list
        self.site = glyc[0].replace('/', '\\')
        self.glycan_list = glycan_list

    def run(self):
        group_li = [1, 2]
        dict_expression = {}
        site_Both_Group = {}
        site_final = {}

        # Volcano Expression Plot 放在外面减少程序资源负担
        # Volcano Expression Plot
        # 指定文件名（硬编码）
        expression_level_file = "output_Volcano.xlsx"  # 在这里指定文件名
        # 构建完整的文件路径
        file_path = pathlib.Path(self.expression_level, expression_level_file)
        # 检查文件是否存在
        if os.path.exists(file_path):
            final_signal.write(f"File '{expression_level_file}' is in '{self.expression_level}'.")
        else:
            final_signal.write(f"File '{expression_level_file}' is not in '{self.expression_level}'.")
        df = pd.read_excel(file_path, engine='openpyxl', sheet_name='Differential Expression')
        # 获取列名（表头）
        columns = list(df.columns)  # 将表头转换为列表
        for pr in self.protein_list:
            # 逐行生成矩阵，每个矩阵包含表头和对应行数据
            for index, row in df.iterrows():
                # 将表头和对应行内容合并为一个矩阵
                matrix = np.array([columns, row.values])
                if row['ID'] == pr:
                    dict_expression[pr] = matrix
            # Site
            tem = {}
            for group in group_li:
                Site_Route1 = pathlib.Path(self.site, f'Figure_GlycanSites_Group{group}')
                file_pattern = re.compile(rf"{group}_(\d+)_Glycan_sites_{pr}\.xlsx")
                xlsx_files = [f for f in os.listdir(Site_Route1) if f.endswith('.xlsx')]
                matching_files = [(int(match.group(1)), file_name) for file_name in xlsx_files if
                                  (match := file_pattern.match(file_name))]
                matching_files.sort(key=lambda x: x[0])
                Site_Total = [site[1] for site in matching_files]
                tem[group] = Site_Total
            site_Both_Group[pr] = tem
        vector_p = self.glycan_vector(self.glycan_list)
        glycan_final = vector_p
        for Protein, Site_group in site_Both_Group.items():
            combined_data_site = self.site_vector(Site_group, Protein)
            site_final[Protein] = combined_data_site
            # dict_expression = vector_e
        return dict_expression, glycan_final, site_final

    def js_divergence(self, p, q):
        """计算两个概率分布的Jensen-Shannon散度"""
        # 处理零输入
        if np.all(p == 0) or np.all(q == 0):
            return np.nan

        # 添加平滑处理
        epsilon = 1e-10
        _p = (p + epsilon) / (np.sum(p) + epsilon * len(p))
        _q = (q + epsilon) / (np.sum(q) + epsilon * len(q))

        _m = 0.5 * (_p + _q)
        try:
            return 0.5 * (distance.entropy(_p, _m) + distance.entropy(_q, _m))
        except:
            return np.nan

    def process_matrix_stats(self, matrix, prefix):
        """增强版矩阵统计量聚合"""
        stat_index = {}
        if matrix is None or len(matrix) == 0:
            # 初始化所有特征为NaN
            stat_index.update({f'{prefix}_OD_js': np.nan, f'{prefix}_OD_range_ratio': np.nan})
            # return stats for test

        # ================== Observed Difference处理 ==================
        if len(matrix) != 0:
            od = []
            for site_glycan, item0 in matrix.items():
                observed_diff = item0['Observed Difference']  # None
                if observed_diff is not None and not np.isnan(observed_diff) and observed_diff != 0:
                    od.append(item0['Observed Difference'])
            od = np.array(od)
            if len(od) != 0:
                # 基础统计量
                stat_index.update({
                    f'{prefix}_OD_mean': od.mean(),
                    f'{prefix}_OD_median': np.median(od),
                    f'{prefix}_OD_std': od.std(),
                    f'{prefix}_OD_range': od.max() - od.min(),
                    f'{prefix}_OD_iqr': np.percentile(od, 75) - np.percentile(od, 25),  # 四分位距
                    f'{prefix}_OD_cv': (od.std() / od.mean() if od.mean() != 0 else np.nan)  # 变异系数
                })

                # 极差与中位数比值（处理除零）
                median_val = np.median(od)
                range_ratio = (od.max() - od.min()) / median_val if median_val != 0 else np.nan
                stat_index[f'{prefix}_OD_range_ratio'] = range_ratio

                # Jensen-Shannon散度
                try:
                    hist_real, bins = np.histogram(od, bins=10, density=True)
                    bin_width = bins[1] - bins[0]
                    hist_real = hist_real * bin_width
                    median_val = np.median(od)
                    ref_dist = stats.norm.pdf(bins[:-1], loc=median_val, scale=od.std())
                    ref_dist = ref_dist / ref_dist.sum()
                    stat_index[f'{prefix}_OD_js'] = self.js_divergence(hist_real, ref_dist)
                except Exception as e:
                    print(f'Error: {e}, return Jensen-Shannon divergence as nan')
                    stat_index[f'{prefix}_OD_js'] = np.nan
            else:
                stat_index.update({
                    f'{prefix}_OD_iqr': np.nan,
                    f'{prefix}_OD_cv': np.nan,
                    f'{prefix}_OD_range_ratio': np.nan,
                    f'{prefix}_OD_js': np.nan
                })
        else:
            stat_index.update({
                f'{prefix}_OD_iqr': np.nan,
                f'{prefix}_OD_cv': np.nan,
                f'{prefix}_OD_range_ratio': np.nan,
                f'{prefix}_OD_js': np.nan
            })

        # ================== p-value处理（Fisher合并法） ==================
        if len(matrix) != 0:
            p = []
            for site_glycan, item0 in matrix.items():
                p_value = item0['p-value']
                if p_value is not None and not np.isnan(p_value) and p_value != 0:
                    p.append(p_value)
            p_values = np.array(p)
            valid_p = p_values[p_values > 0]  # 排除0和NaN
            if len(valid_p) > 0:
                # Fisher合并计算
                chi_sq = -2 * np.sum(np.log(valid_p))
                df = 2 * len(valid_p)
                combined_p = chi2.sf(chi_sq, df)
                stat_index.update({
                    f'{prefix}_p_combined': combined_p,
                    f'{prefix}_p_sig_count': (p_values < 0.05).mean()
                })
            else:
                stat_index.update({f'{prefix}_p_combined': np.nan, f'{prefix}_p_sig_count': 0})
        else:
            stat_index.update({f'{prefix}_p_combined': np.nan, f'{prefix}_p_sig_count': 0})

        # ================== 95% CI/HDI处理 ==================
        if len(matrix) != 0:
            CI = []
            for site_glycan, item0 in matrix.items():
                # 使用 get 避免 KeyError，并检查键是否存在
                ci_data = item0.get('95% CI')  # 更安全的取值方式
                if ci_data is not None:
                    # 检查 ci_data 是否为列表或元组
                    if isinstance(ci_data, (list, tuple)):
                        # 使用列表推导式过滤 None 值
                        filtered_ci = [c_i for c_i in ci_data if c_i is not None]
                        # 仅当过滤后的列表非空时添加到 CI
                        if filtered_ci:
                            CI.append(filtered_ci)
            try:
                ci_lower = []
                ci_upper = []
                ci_width = []
                for ci_data in CI:
                    try:
                        lower = ci_data[0]
                        upper = ci_data[1]
                        width = upper - lower
                        ci_lower.append(lower)
                        ci_upper.append(upper)
                        ci_width.append(width)
                    except:
                        pass
                ci_lower = np.array(ci_lower)
                ci_upper = np.array(ci_upper)
                ci_width = np.array(ci_width)
                # 基础特征
                stat_index.update({
                    f'{prefix}_HDI/CI_lower_mean': ci_lower.mean(),
                    f'{prefix}_HDI/CI_upper_mean': ci_upper.mean(),
                    f'{prefix}_HDI/CI_width_mean': ci_width.mean(),
                    f'{prefix}_HDI/CI_width_std': ci_width.std(),
                    f'{prefix}_HDI/CI_non_zero': ((ci_lower > 0) | (ci_upper < 0)).sum()
                })
                stat_index[f'{prefix}_HDI/CI_lower_std'] = ci_lower.std()
                stat_index[f'{prefix}_HDI/CI_upper_std'] = ci_upper.std()
                stat_index[f'{prefix}_HDI/CI_lower_range'] = ci_lower.max() - ci_lower.min()
                stat_index[f'{prefix}_HDI/CI_upper_range'] = ci_upper.max() - ci_upper.min()
            except Exception as e:
                pass  # 保持原有NaN处理
        return stat_index

    # singal glycosylation matrixes combination
    def prepare_data(self, matrices, index=1):
        fc_matrix, glycan_matrix, site_matrix, byonic_matrix = matrices

        fc_matrix_after = safe_flatten(fc_matrix)
        fc_matrix_dic = {}
        if fc_matrix_after is not None and fc_matrix_after.shape == fc_matrix.shape:
            raise ValueError('Error: Data type and structure errors (Fold-change/P-Value matrix)')
        elif fc_matrix_after is not None and fc_matrix_after.shape != fc_matrix.shape:
            keys = fc_matrix_after[:3]
            values = fc_matrix_after[3:]
            fc_matrix_dic = {keys[i]: values[i] for i in range(len(values))}
        # 生成增强特征
        site_stats = self.process_matrix_stats(site_matrix, 'Site')
        glycan_stats = self.process_matrix_stats(glycan_matrix, 'Glycan')
        # 构建最终数据集
        combined_data = pd.DataFrame({
            'Fc': fc_matrix_dic.get('Fold Change-Logarithmic form', np.nan),
            'p-Value': fc_matrix_dic.get('p-Value-Logarithmic form', np.nan),
            **glycan_stats,
            **site_stats,
            'Byonic_bootstrap_OD': byonic_matrix.get('Mean Difference', np.nan),
            'Byonic_bootstrap_CI_Lower': min(byonic_matrix['95% CI']) if '95% CI' in byonic_matrix and len(
                byonic_matrix['95% CI']) > 0 else None,
            'Byonic_bootstrap_CI_Upper': max(byonic_matrix['95% CI']) if '95% CI' in byonic_matrix and len(
                byonic_matrix['95% CI']) > 0 else None
        }, index=[index])
        return combined_data

    def xgboost_analysis(self, data_input, outputpath, ForAnnotation, plot=True):
        # 数据预处理部分保持不变
        fc = data_input['Fc']
        pv = data_input['p-Value']
        li1 = pd.concat([fc, pv], axis=1)
        li1.columns = ['Fold Change-Logarithmic form', 'p-Value-Logarithmic form']

        y_XG_df = pd.DataFrame(li1)
        y_XG_df = y_XG_df.apply(pd.to_numeric, errors='coerce')

        print("y_XG_df before filling NaNs:\n", y_XG_df)

        if y_XG_df.isnull().values.any():
            print("Warning: y_XG contains NaN values. Filling with mean.")
            imputer = SimpleImputer(strategy='mean')
            y_XG_df = pd.DataFrame(imputer.fit_transform(y_XG_df), columns=y_XG_df.columns)

        print("y_XG_df after filling NaNs:\n", y_XG_df)

        y_XG_df['Product'] = y_XG_df['Fold Change-Logarithmic form'] * y_XG_df['p-Value-Logarithmic form']
        y_XG_df['Sym_Log_Product'] = y_XG_df['Product'].apply(sym_log_transform)
        y_XG = y_XG_df['Sym_Log_Product'].to_numpy()

        if np.isnan(y_XG).any():
            raise ValueError("y_XG contains NaN!")

        # 特征处理（移除标准化）
        X_XG = data_input.drop(columns=['Fc', 'p-Value', 'ID'])
        X_XG = X_XG.apply(pd.to_numeric, errors='coerce')

        if X_XG.isnull().values.any():
            print("Warning: X_XG contains NaN values. Filling with mean.")
        imputer = SimpleImputer(strategy='mean')
        X_imp = imputer.fit_transform(X_XG)  # 直接使用填充后的数据，不标准化
        X_XG_processed = pd.DataFrame(X_imp, columns=X_XG.columns)

        # 保存预处理对象（移除scaler）
        self.preprocessor = {
            'imputer': imputer,
            'feature_names': X_XG.columns.tolist()
        }

        # XGBoost参数网格
        param_grid = {
            'n_estimators': [100, 200],
            'learning_rate': [0.05, 0.1],
            'max_depth': [3, 5],
            'subsample': [0.8, 1.0]
        }

        # 模型训练
        model = GridSearchCV(
            estimator=XGBRegressor(random_state=42, n_jobs=-1),
            param_grid=param_grid,
            cv=5,
            scoring='neg_mean_squared_error',
            verbose=1
        ).fit(X_XG_processed, y_XG)

        best_params = model.best_params_
        best_mse = -model.best_score_
        print(f"Best Params: {best_params}, Cross-validated MSE: {best_mse:.4f}")

        # 特征重要性
        coefs = pd.Series(
            model.best_estimator_.feature_importances_,
            index=self.preprocessor['feature_names']
        ).sort_values(ascending=False)

        # 蛋白质影响力评分
        protein_impact = {}
        try:
            for protein in data_input['ID'].unique():
                protein_data = data_input[data_input['ID'] == protein]
                X_protein = protein_data.drop(columns=['Fc', 'p-Value', 'ID'])
                X_imp_protein = self.preprocessor['imputer'].transform(X_protein)
                score = model.best_estimator_.predict(X_imp_protein).mean()
                protein_impact[protein] = score
        except Exception as e:
            raise RuntimeError(f"Error in protein impact calculation: {str(e)}")

        # 保存结果
        impact_df = pd.DataFrame({
            'Protein': list(protein_impact.keys()),
            'Impact_Score': list(protein_impact.values())
        }).sort_values('Impact_Score', ascending=False)

        impact_path = os.path.join(outputpath, 'protein_impact_scores.xlsx')
        impact_df.to_excel(impact_path, index=False)
        print(f"Impact Score saved at: {impact_path}")

        # 可视化
        if plot:
            plt.figure(figsize=(12, 8))
            sns.barplot(x='Impact_Score', y='Protein', data=impact_df.head(20), palette='viridis')
            plt.title('Top 20 Proteins by Feature Impact')
            plt.tight_layout()
            plt.savefig(os.path.join(outputpath, 'key_proteins_impact.png'), dpi=1000)
            plt.close()


        # 模型评估
        y_pred = model.best_estimator_.predict(X_XG_processed)
        mse = mean_squared_error(y_XG, y_pred)
        rmse = np.sqrt(mse)
        r2 = r2_score(y_XG, y_pred)

        print('\n=== Model Evaluation ===')
        print(f"MSE: {mse:.4f}")
        print(f"RMSE: {rmse:.4f}")
        print(f"R²: {r2:.4f}")

        metrics_df = pd.DataFrame({
            'Metric': ['MSE', 'RMSE', 'R²'],
            'Value': [mse, rmse, r2]
        })
        metrics_path = os.path.join(outputpath, 'XGBoost_metrics.xlsx')
        metrics_df.to_excel(metrics_path, index=False)

        # 二级Annotaion
        ori = pathlib.Path(outputpath)
        annotation_dir = ori / 'Annotation'
        annotation_dir.mkdir(parents=True, exist_ok=True)
        Annotation.main(impact_path=impact_path,path_list=ForAnnotation,path_save_plot=annotation_dir)
        return coefs, impact_df, impact_path

    def format_glycan_structure(self, glycan_structure_list):
        formatted_glycans = []
        sugar_code_map = {
            'H': 'Hex',
            'N': 'HexNAc',
            'F': 'Fuc',
            'A': 'NeuAc',
            'G': 'NeuGc',
        }
        sugar_order = ['Hex', 'HexNAc', 'Fuc', 'NeuAc', 'NeuGc']
        if type(glycan_structure_list) == list:
            pass
        elif type(glycan_structure_list) == str:
            glycan_structure_list = [glycan_structure_list]
        for item in glycan_structure_list:
            if isinstance(item, list):
                glycan_str = item[0] if item else ''
            else:
                glycan_str = str(item)

            glycan_dict = {}

            # Type: {Hex:3; Fuc:1}
            if glycan_str.startswith('{'):
                matches = re.findall(r'([A-Za-z0-9]+):(\d+)', glycan_str)
                for sugar, count in matches:
                    # 统一命名：Neu5Ac/NeuAc → NeuAc；Neu5Gc/NeuGc → NeuGc
                    sugar = sugar.replace('Neu5Ac', 'NeuAc').replace('Neu5Gc', 'NeuGc')
                    glycan_dict[sugar] = glycan_dict.get(sugar, 0) + int(count)

            # Type: H1N1A2
            else:
                matches = re.findall(r'([A-Z])(\d+)', glycan_str)
                for code, count in matches:
                    sugar = sugar_code_map.get(code)
                    if sugar:
                        glycan_dict[sugar] = glycan_dict.get(sugar, 0) + int(count)

            components = []
            for sugar in sugar_order:
                if sugar in glycan_dict:
                    components.append(f"{sugar}({glycan_dict[sugar]})")
            formatted_glycans.append(''.join(components))

        return formatted_glycans

    # Bootstrap 显著性分析函数
    def bootstrap_significance(self, entries, alpha=0.05):
        """Run bootstrap significance, accepting either a DataFrame with 'Intensity' & 'Group' or
        a dict mapping group -> list of intensity values."""
        # Convert dict input to DataFrame
        if isinstance(entries, dict):
            records = []
            for grp, vals in entries.items():
                # Ensure list of values
                if not isinstance(vals, (list, np.ndarray)):
                    vals = [vals] if vals is not None else []
                for v in vals:
                    records.append({'Intensity': v, 'Group': grp})
            df = pd.DataFrame(records)
        elif isinstance(entries, pd.DataFrame):
            df = entries.copy()
            # Ensure required columns exist
            if 'Intensity' not in df.columns or 'Group' not in df.columns:
                raise ValueError("DataFrame must contain 'Intensity' and 'Group' columns.")
        else:
            raise ValueError("Entries must be a DataFrame or a dict of group->values.")

        # Require at least two groups and samples
        if df['Group'].nunique() < 2 or df.shape[0] < 2:
            return None, np.nan, np.nan, np.nan

        grp_sizes = df.groupby('Group')['Intensity'].count()
        grp_means = df.groupby('Group')['Intensity'].mean()
        grp_vars = df.groupby('Group')['Intensity'].var()

        observed_diff, _ = compute_effect_size(grp_means, grp_vars, grp_sizes)
        final_signal.write('Deploying bootstrap/bayesian iteration tasks on Ray cluster....')
        worker_ray_path = pathlib.Path(resource_path("resources")) / "result_cache_2.pickle"
        if worker_ray_path.exists():
            worker_ray_path.unlink()

        connection_docker_2(entries)
        container = docker_ray_manager_2.start_ray_container()

        # log = pathlib.Path(current_directory, "docker_ray",'task1_complete.txt')
        while not worker_ray_path.exists():
            time.sleep(3)
        final_signal.write('Bootstrap tasks of glycan or glycosylation site completed')

        with open(worker_ray_path, 'rb') as handle:
            bootstrap_diffs = pickle.load(handle)

        final_signal.write(f'Finished the iteration tasks on Ray cluster in docker. Output: {bootstrap_diffs}')
        docker_ray_manager_2.stop_ray_container(container)
        if bootstrap_diffs and bootstrap_diffs[0][0] == 'bootstrap':
            # 1. 初始化
            bootstrap_effects = []
            perm_effects = []

            # 2. 遍历所有迭代结果
            for item in bootstrap_diffs:
                tag, eff, *_, perm_eff = item
                # 只处理标记为 'bootstrap' 的
                if tag != 'bootstrap':
                    continue

                # 收集 eff
                if not np.isnan(eff):
                    bootstrap_effects.append(eff)
                # 收集 perm_eff
                if not np.isnan(perm_eff):
                    perm_effects.append(perm_eff)

            # 3. 如果没有有效的 bootstrap 结果，则直接返回
            if len(bootstrap_effects) == 0:
                return observed_diff, np.nan, (np.nan, np.nan), np.nan

            # 4. 计算 p-value 与 BCa 置信区间
            p_value, (lower, upper) = calculate_pvalue_and_ci(
                observed_diff,
                np.array(bootstrap_effects),
                df,
                alpha
            )

            # 5. 置换检验 p 值
            perm_p = (np.sum(np.abs(perm_effects) >= np.abs(observed_diff)) + 1) \
                     / (len(perm_effects) + 1)

            # 6. 打日志
            final_signal.write(f"Bootstrap p: {p_value:.4f}, Permutation p: {perm_p:.4f}\n")
            final_signal.write(f"95% BCa CI: [{lower:.3f}, {upper:.3f}]\n")

            return observed_diff, p_value, (lower, upper), perm_p
        else:
            bayesian_effect = bootstrap_diffs
            # —— Bayesian 分支 ——
            # 从 bayesian_effect 获取：mean_diff, (ci_low, ci_high), bayes_p
            tag, dif, low, high, bayes_p, _ = bayesian_effect[0]
            assert tag == 'bayesian', "Unexpected tag in Bayesian branch"
            # 将 perm_p 固定为 np.nan，以便下游不混淆
            final_signal.write(f"Bayesian HDI: [{low:.3f}, {high:.3f}], Bayes p: {bayes_p:.4f}")
            return dif, bayes_p, (low, high), np.nan
        # return observed_diff, p_value, (lower, upper), perm_p

    def site_vector(self, Site_group, Protein):
        # Build DataFrame
        records = []
        for grp, files in Site_group.items():
            for fname in files:
                df = pd.read_excel(pathlib.Path(self.site) / f'Figure_GlycanSites_Group{grp}' / fname,
                                   engine='openpyxl')
                for _, row in df.iterrows():
                    if pd.notna(row.get('Sites')) and pd.notna(row.get('Intens.')):
                        records.append(
                            {'Protein': Protein, 'Site': row['Sites'], 'Intensity': row['Intens.'], 'Group': grp})
        df_sites = pd.DataFrame(records).dropna()
        if df_sites.empty:
            return None

        # For each unique Site perform bootstrap across groups
        results = {}
        for site in df_sites['Site'].unique():
            sub = df_sites[df_sites['Site'] == site]
            if sub['Group'].nunique() < 2 or sub.shape[0] < 2:
                results[site] = {'Observed Difference': None, 'p-value': None, '95% CI': (None, None)}
                continue
            od, pval, ci, _ = self.bootstrap_significance(sub)
            results[site] = {'Observed Difference': od, 'p-value': pval, '95% CI': ci}
        return results

    def glycan_vector(self, glycan_dict):
        """
        Compute bootstrap stats per protein per glycan.
        Input: glycan_dict = { group: [struct_map, intensity_map], ... }
        Returns: { protein_id: { glycan_str: {'Observed Difference', 'p-value', '95% CI'}, ... }, ... }
        """
        results = {}
        # Iterate each group and build per-protein records
        for grp, maps in glycan_dict.items():
            struct_map, int_map = maps
            for file_idx, prot_map in struct_map.items():
                signals = int_map.get(file_idx, {})
                for prot_id, struct in prot_map.items():
                    formatted = self.format_glycan_structure(struct)
                    vals = signals.get(prot_id)
                    vals_list = vals if isinstance(vals, list) else ([vals] if vals is not None else [])
                    if len(formatted) != len(vals_list):
                        final_signal.write(f"Structure/intensity length mismatch for {prot_id} in file {file_idx}")
                        continue
                    # prepare per protein grouping
                    if prot_id not in results:
                        results[prot_id] = {}
                    # accumulate intensities per glycan
                    for gly_str, inten in zip(formatted, vals_list):
                        results[prot_id].setdefault(gly_str, []).append((inten, grp))
        # Now compute bootstrap per protein per glycan
        final = {}
        for prot_id, gly_map in results.items():
            final[prot_id] = {}
            for gly_str, samples in gly_map.items():
                # samples: list of (intensity, group)
                df = pd.DataFrame(samples, columns=['Intensity', 'Group'])
                od, p, ci, _ = self.bootstrap_significance(df)
                final[prot_id][gly_str] = {'Observed Difference': od, 'p-value': p, '95% CI': ci}
        return final


class EnhancedML:
    def __init__(self, data, output_path):
        self.data = data
        self.output_path = output_path
        self.best_model = None
        self.feature_importances = None
        self.pca_results = None

        # 执行完整分析流程
        self.run_analysis()

    def generate_pca_results(self, data):
        """生成包含蛋白质ID的PCA结果DataFrame"""
        import pandas as pd
        # 假设 data 中包含原始ID列
        pca_df = pd.DataFrame(
            self.pca_result,  # 来自 train_optimized_rf 中计算的 self.pca_result
            columns=[f'PC{i + 1}' for i in range(self.pca_result.shape[1])]
        )
        pca_df['Protein'] = data['ID'].values
        return pca_df

    def run_analysis(self):
        """执行完整分析流程"""
        try:
            # 数据预处理
            preprocessed_data = self.preprocess_data()

            # 参数优化与模型训练
            self.train_optimized_rf(preprocessed_data)

            # 特征加权PCA分析（直接使用已计算的PCA结果）
            self.pca_results = self.generate_pca_results(preprocessed_data)

            # 获取蛋白质描述信息
            Pr = list(self.pca_results['Protein'])
            c = []
            for ID in Pr:
                try:
                    o = self.get_info_biopython(ID)
                    comment = o.get("Comments", "No description available")
                except Exception as e:
                    comment = f"Error fetching description: {str(e)}"
                c.append(str(comment))

            # 插入描述列
            self.pca_results.insert(1, 'Description', c)

            # 结果保存与可视化
            self.save_results()
            self.visualize_clustering()  # 确保方法名统一

        except Exception as e:
            print(f"分析流程失败: {str(e)}")
            raise

    def preprocess_data(self):
        """数据预处理"""
        # 清除缺失值超过60%的列
        # 检查Data数据类型
        cleaned_data = self.data.dropna(thresh=int(0.6 * len(self.data)), axis=1)
        # 填充剩余缺失值
        numeric_cols = cleaned_data.select_dtypes(include=np.number).columns
        cleaned_data[numeric_cols] = cleaned_data[numeric_cols].fillna(cleaned_data[numeric_cols].mean())

        return cleaned_data

    def train_optimized_rf(self, data):
        """自动选择聚类数，生成伪标签并分析特征重要性"""
        X = data.drop(columns=['ID'])

        # --- 1. 自动选择最佳聚类数 ---
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score

        k_range = range(2, 11)
        self.best_k = 2  # 保存到 self 属性
        best_score = -1

        for k in k_range:
            kmeans = KMeans(n_clusters=k, random_state=42)
            labels = kmeans.fit_predict(X)

            if len(np.unique(labels)) < 2:
                continue

            score = silhouette_score(X, labels)
            print(f"K={k} | Silhouette Score={score:.3f}")

            if score > best_score:
                best_score = score
                self.best_k = k  # 更新 self.best_k
                self.final_kmeans = kmeans  # 保存最佳模型到 self

        print(f"自动选择聚类数: K={self.best_k}")

        # --- 2. 生成伪标签 ---
        pseudo_labels = self.final_kmeans.fit_predict(X)

        # --- 3. 随机森林特征重要性 ---
        from sklearn.ensemble import RandomForestClassifier
        rf = RandomForestClassifier(random_state=42)
        rf.fit(X, pseudo_labels)

        self.feature_importances = pd.Series(
            rf.feature_importances_,
            index=X.columns
        ).sort_values(ascending=False)

        # --- 4. 动态权重调整PCA ---
        # 计算累计重要性，确定保留维度
        cumulative_importance = self.feature_importances.cumsum()
        n_components = np.argmax(cumulative_importance >= 0.8) + 1  # 保留80%重要性的特征

        # 确保至少保留2个维度
        n_components = max(2, n_components) if n_components > 0 else 2

        # 特征加权
        weighted_X = X * self.feature_importances.values

        from sklearn.decomposition import PCA
        self.pca = PCA(n_components=n_components)
        self.pca_result = self.pca.fit_transform(weighted_X)

        print("Importance ranking of features:\n", self.feature_importances)
        print(f"PCA retained dimensions: {n_components}，PCA variance explained: {self.pca.explained_variance_ratio_}")

    def visualize_clustering(self):
        """可视化聚类结果与PCA投影"""
        import matplotlib.pyplot as plt

        # 确保已定义关键变量
        if not hasattr(self, 'final_kmeans') or not hasattr(self, 'best_k'):
            raise ValueError("请先运行 train_optimized_rf() 生成聚类结果")

        plt.figure(figsize=(12, 5))

        # 子图1：PCA投影
        plt.subplot(121)
        scatter = plt.scatter(
            self.pca_result[:, 0],
            self.pca_result[:, 1],
            c=self.final_kmeans.labels_,
            cmap='viridis',
            edgecolor='k'
        )
        plt.xlabel(f'PC1 ({self.pca.explained_variance_ratio_[0]:.1%})')
        plt.ylabel(f'PC2 ({self.pca.explained_variance_ratio_[1]:.1%})')
        plt.title(f'K={self.best_k} 的PCA投影')

        # 子图2：特征重要性
        plt.subplot(122)
        self.feature_importances.plot(kind='barh')
        plt.title('Importance ranking of features')
        plt.tight_layout()
        plt.show()

    def save_results(self):
        """保存关键结果"""
        # 保存特征重要性
        self.feature_importances.to_csv(
            os.path.join(self.output_path, 'feature_importances.csv'))

        # 保存分类结果
        self.pca_results.to_excel(
            os.path.join(self.output_path, 'protein_classification.xlsx'),
            index=False)

        # 保存模型参数
        # with open(os.path.join(self.output_path, 'best_params_RFClassify.pkl'), 'wb') as f:
        #     pickle.dump(self.best_model.get_params(), f)

    def get_info_biopython(self, uniprot_id):
        try:
            handle = ExPASy.get_sprot_raw(uniprot_id)
            record = SwissProt.read(handle)
            comments = record.comments
            # 处理 comments 结构
            if comments:
                first_comment = comments[0]
                if isinstance(first_comment, list):
                    first_comment = "; ".join(map(str, first_comment))
            else:
                first_comment = "No comment"
            return {"Comments": first_comment}
        except Exception as e:
            return {"Comments": f"Error: {str(e)}"}


def hide_console_window():
    """隐藏当前进程的控制台窗口（仅限 Windows）。"""
    if os.name == "nt":  # 仅限 Windows 环境
        kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        user32 = ctypes.WinDLL('user32', use_last_error=True)

        SW_HIDE = 0
        hwnd = kernel32.GetConsoleWindow()
        if hwnd:
            user32.ShowWindow(hwnd, SW_HIDE)


def safe_flatten(matrix):
    if isinstance(matrix, pd.DataFrame):
        data_cleaned = matrix.drop(columns=['Unnamed: 0'])
        return data_cleaned.values.flatten()
    elif isinstance(matrix, np.ndarray):
        data_cleaned = matrix[:, 1:]  # Remove the first column
        return data_cleaned.flatten()
    else:
        return matrix


# 定义对称化对数变换函数
def sym_log_transform(x):
    return np.sign(x) * np.log1p(np.abs(x))


class PCA_Glycosylation_Analysis:
    def __init__(self, byonic1, byonic2, Total, T, output, ForAnnotation):

        # glycoprotein_data = {
        #     'Glycoprotein1': [fc_matrix_1, glycan_matrix_1, site_matrix_1, byonic_matrix_1],
        #     'Glycoprotein2': [fc_matrix_2, glycan_matrix_2, site_matrix_2, byonic_matrix_2],
        # }
        # cfg.Engine.put("Ray")
        # start_ray_cluster("rayproject/ray:latest-py312")
        # wait_for_ray_cluster()
        byonic = byonic_r(byonic1, byonic2)
        protein_list = byonic.protein_list
        glycan_list = []
        result_by = byonic.bootstrap_analysis()
        try:
            file_to_find1 = "GlycReSoft_Tem_file_for_FinalAnalysis.pkl"
            file_to_find2 = "MetaMorpheus_Tem_file_for_FinalAnalysis.pkl"
            for file_name in os.listdir(T[0].replace('/', '\\')):
                file_path = pathlib.Path(T[0].replace('/', '\\'), file_name)
                if os.path.isfile(file_path) and file_name == file_to_find1:
                    final_signal.write(f'{file_path} is found')
                    with open(file_path, 'rb') as f:
                        glycan_list = pickle.load(f)
                elif os.path.isfile(file_path) and file_name == file_to_find2:
                    final_signal.write(f'{file_path} is found')
                    with open(file_path, 'rb') as f:
                        glycan_list = pickle.load(f)
        except FileNotFoundError as e:
            final_signal.write(e)
        if len(glycan_list) != 0:
            VectorMake = Vector(T, Total, protein_list, glycan_list)
        else:
            raise FileNotFoundError
        glycoprotein_data = {}
        # every matrix from different glycoproteins except byonic
        fc_matrix, glycan_matrix, site_matrix = VectorMake.run()

        for pr in protein_list:
            try:
                a = fc_matrix[pr]
            except Exception as e:
                final_signal.write(f'For {pr},Fc matrix has occured issue, {e}')
                a = None
            try:
                b = glycan_matrix[pr]
            except Exception as e:
                final_signal.write(f'For {pr},glycan matrix has occured issue, {e}')
                b = None
            try:
                c = site_matrix[pr]
            except Exception as e:
                final_signal.write(f'For {pr},site matrix has occured issue {e}')
                c = None
            try:
                d = result_by[pr]
            except Exception as e:
                final_signal.write(f'For {pr},glycan matrix has occured issue, {e}')
                d = None

            try:
                glycoprotein_data[pr] = [a, b, c, d]
            except Exception as e:
                final_signal.write(e)
                glycoprotein_data[pr] = None
                continue

        # 将glycoprotein_data转为Ridge分析所需的DataFrame格式
        combined_data = pd.DataFrame()

        for pr in protein_list:
            # 提取每个糖蛋白的数据并格式化为DataFrame
            data_for_ridge = VectorMake.prepare_data(glycoprotein_data[pr])
            data_for_ridge['ID'] = pr
            # 将每个蛋白质的数据合并到总的DataFrame中
            combined_data = pd.concat([combined_data, data_for_ridge], ignore_index=True)

        # 调用ridge_analysis并传入格式化后的DataFrame
        df_fillNan = combined_data.fillna(0)
        coef, impact_df,impact_path = VectorMake.xgboost_analysis(df_fillNan, output, ForAnnotation)
        # Ensure 'coef' is converted to a dictionary if it's a pandas Series
        # coef = coef.to_dict()
        # # Convert dictionary to DataFrame
        # df = pd.DataFrame(list(coef.items()), columns=['Feature', 'Value'])  # Each key-value pair becomes a row
        # # Write DataFrame to Excel file
        # output_file = 'glycoproteins_data_coef.xlsx'
        # path = os.path.join(output, output_file)
        # df.to_excel(path, index=False)
        print('--------------------------')
        print('Coefficient of XGBoost Analysis')
        print(f"{coef}")
        print('--------------------------')
        print(f'{impact_df}')

        tem = pathlib.Path(resource_path('resources'))
        # os.remove(tem / 'task2_complete.txt')
        os.remove(tem / 'task1_complete.txt')

        # 调用RandomForest并传入格式化后的字典
        EnhancedML(df_fillNan, output)
        final_signal.write('Finished all jobs')