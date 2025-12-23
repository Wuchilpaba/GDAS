import traceback
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
import pandas
from math import log10,log2
import math
from statsmodels.stats.multitest import fdrcorrection

matplotlib.use('agg') # 后端渲染或者用'svg'
mpl.rcParams['font.family'] = 'Times New Roman'
class Statistic:
    def Volcano(self):
        dic1 = self.dic1
        dic2 = self.dic2
        excel_file_path = f'{self.save}\\output_Volcano.xlsx'
        lis_indx = sorted(list(dic1.keys()))
        Fc = []
        p = []
        for item in lis_indx:
            Fc.append(log2(dic1[item]))
            p.append((-1)*log10(dic2[item]))
        Frame = pandas.DataFrame({'ID':lis_indx,'Fold Change-Logarithmic form':Fc, 'p-Value-Logarithmic form':p})

        Frame.to_excel(excel_file_path,sheet_name='Differential Expression')

        # 设置阈值,FC=1.5,P=0.05
        threshold_x = log2(self.FcValue) # log2(1.5)
        threshold_x_neg = -log2(self.FcValue) # -log2(1.5)
        threshold_y = -log10(self.PValue) # -log10(0.05)
        # x = FC dic1
        # y = P-Value dic2
        try:
            x = Fc
            y = p
            # 根据条件将数据点分组，并使用不同颜色绘制
            colors = []
            for xi, yi in zip(x, y):
                if xi > threshold_x and yi > threshold_y:
                    colors.append('red')
                elif xi < threshold_x_neg and yi > threshold_y:
                    colors.append('blue')
                else:
                    colors.append('gray')
            print(list(zip(x, y)))

            # 统计x轴大于Fc，y轴大于p的数据数量和x轴小于-p，y轴大于2的数据数量
            count_UP = len([1 for xi, yi in zip(x, y) if xi > threshold_x and yi > threshold_y])
            count_DOWN = len([1 for xi, yi in zip(x, y) if xi < threshold_x_neg and yi > threshold_y])

            # 绘制火山图
            plt.figure(figsize=(8, 6))
            plt.scatter(x, y, color=colors, alpha=0.8, s=35)  # 使用散点图绘制
            max_abs_x = max(abs(xi) for xi in x)  # 计算 x 列表中最大绝对值
            plt.xlim(-math.ceil(max_abs_x) - 1, math.ceil(max_abs_x) + 1)
            # plt.xlim(-5,5)
            # plt.ylim(-2,10)
            # 添加标注信息
            plt.text(0.85, 0.9, f'UP: {count_UP}', color='red', fontsize=10, transform=plt.gcf().transFigure)
            plt.text(0.13, 0.9, f'DOWN: {count_DOWN}', color='blue', fontsize=10, transform=plt.gcf().transFigure)

            # 添加虚线表示阈值
            plt.axvline(x=threshold_x, linestyle='--', color='gray')  # 添加第一个垂直虚线
            plt.axvline(x=threshold_x_neg, linestyle='--', color='gray')  # 添加第二个垂直虚线
            plt.axhline(y=threshold_y, linestyle='--', color='gray')  # 添加水平虚线

            plt.xlabel('log2(Fold Change)')  # 设置横坐标标签
            plt.ylabel('-log10(P-Value)')  # 设置纵坐标标签
            plt.title('Differential Expression Analysis')  # 设置标题
            plt.grid(True, color='lightgrey')  # 显示网格线并设置颜色为浅灰色
            plt.grid(False)  # 不显示网格线
            # plt.autoscale()
            plt.savefig(f"{self.save}\\Vocanlo.png", dpi=1000)

        except Exception as e:
            traceback.print_exc()
            print(e)
        finally:
            print('Differential Expression Analysis Plotting is Finished.')
            return Fc,p,lis_indx
    def __init__(self, dic1, dic2, save, FcValue, PValue):
        p_ID = {}
        FC_ID = {}
        super().__init__()
        self.save = save
        self.FcValue = FcValue
        self.PValue = PValue
        try:
            for key,value in dic1.items():
                lst1 = value
                lst2 = dic2[key]

                s1 = 0
                for s in lst1:
                    s1 += s
                s2 = 0
                for s in lst2:
                    s2 += s
                Fc = 1
                try:
                    a1 = s1/len(lst1)
                    a2 = s2/len(lst2)
                    Fc = a1 / a2
                except ZeroDivisionError:
                    print(f'ZeroDivisionError with {key}')
                    pass

                FC_ID[key] = Fc
                # 独立的两个样本 t 检验
                sample1 = np.asarray(lst1)  # 转换为 numpy 数组
                sample2 = np.asarray(lst2)

                # 检查方差齐性
                levene_test = stats.levene(sample1, sample2)
                print(levene_test)  # 打印 Levene's Test 结果

                # 方差齐性的检验结果
                if levene_test.pvalue > 0.05:
                    # 方差齐性时进行的 t 检验
                    t_test_result = stats.ttest_ind(sample1, sample2)
                else:
                    # 方差不齐时进行的 t 检验
                    t_test_result = stats.ttest_ind(sample1, sample2, equal_var=False)
                sta = t_test_result.statistic
                p = t_test_result.pvalue
                # 打印 t 检验结果
                print(f"{key}_statistic:", sta)
                print(f"{key}_pvalue:", p)
                p_ID[key] = p
                # 基于 p 值判断原假设的接受或拒绝
                if t_test_result.pvalue < 0.05:
                    print(f"{key}: Reject the null hypothesis. There is a significant correlation between the two datasets.")
                else:
                    print(f"{key}: Fail to reject the null hypothesis. There is no significant correlation between the two datasets.")
        except Exception as e:
            print(e)
            traceback.print_exc()
        self.dic1 = FC_ID
        self.dic2 = p_ID