# Mode: After GlycReSoft
import matplotlib
matplotlib.use('agg') # 后端渲染或者用'svg'
import os,logging,openpyxl,sys
import re
import threading
import pandas as pd
from AnalysisSignalModule import AnalysisSignal
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log,log2
from collections import defaultdict
# import resource_path

analysis_signal = AnalysisSignal()
class Analysis_GlycReSoft(threading.Thread):
    def stop(self):
        """设置停止标志为True，用于停止分析进程。"""
        self._stop_event.set()

    def stopped(self):
        """检查是否有停止请求。"""
        return self._stop_event.is_set()

    def check_stop(self):
        if self.stopped():
            print("Analysis was stopped!")
            return True
        return False
    @staticmethod
    def find_index_with_regex(string_list, pattern):
        """
        在字符串列表中使用正则表达式查找匹配的元素，并返回其索引。

        Args:
        - string_list (list): 包含字符串的列表。
        - pattern (str): 要匹配的正则表达式模式。

        Returns:
        - int or None: 匹配的元素的索引，如果未找到匹配项则返回None。
        """
        for index, string in enumerate(string_list):
            if re.search(pattern, string):
                return index
        return None  # 如果没有找到匹配项，则返回None。
    @staticmethod
    def find_duplicates(lst):
        element_dict = {}  # 用于存储元素及其索引的字典
        duplicates = {}  # 存储重复元素及其所有索引的字典
        unique_elements_dict = {}  # 存储不重复元素及其索引的字典
        for i, item in enumerate(lst):
            # 如果元素已经在字典中，说明是重复元素
            if item in element_dict:
                if item not in duplicates:
                    duplicates[item] = [element_dict[item]]
                # 检查当前索引是否已经在列表中
                if i not in duplicates[item]:
                    duplicates[item].append(i)
                    # 从不重复字典中删除已经出现在重复字典中的元素
                    if item in unique_elements_dict:
                        del unique_elements_dict[item]
            else:
                # 将元素及其索引添加到字典中
                element_dict[item] = i
                # 如果元素已经在不重复字典中，不再更新其索引
                if item not in unique_elements_dict:
                    unique_elements_dict[item] = i
        return duplicates, unique_elements_dict
    def __init__(self, file_psm1, file_psm2, Save_path, fasta_input, Multiple_selection):
        super().__init__()
        self.path = {1:file_psm1, 2:file_psm2}
        self.Save = Save_path
        self.fo = []
        self.Both_Peptide = {}
        self.both_intensity_all = {}
        self.Both_ID_list = []
        self.fasta_input = fasta_input
        self.Multiple_selection = Multiple_selection
        self.Ratio3 = {}
        self._stop_event = threading.Event()

        analysis_signal.write(f"Starting analysis\nCSV Group1: {file_psm1}, \nCSV Group2: {file_psm2}, \nSave path: {Save_path}")
        logging.info(f"Starting analysis with \nCSV Group1: {file_psm1}, \nCSV Group2: {file_psm2}, \nSave path: {Save_path}")


    def run(self):
        try:
            self.final = {}
            for indx,item in self.path.items():
                self.final[indx] = {}
                self.Both_Peptide[indx] = {}
                self.both_intensity_all[indx] = {}
                sites_dic = {}
                sites_signal_dic = {}
                glycan_dic = {}
                glycan_signal_dic = {}
                glucose = {}
                glucose_use_signal = {}
                for ifl,file in enumerate(item):
                    f = str(file)
                    form = pd.read_csv(f, encoding='utf-8')
                    logging.info(form)
                    analysis_signal.write(form)
                    glucose[ifl] = {}
                    glucose_use_signal[ifl] = {}
                    sites_dic[ifl] = {}
                    sites_signal_dic[ifl] = {}
                    glycan_dic[ifl] = {}
                    glycan_signal_dic[ifl] = {}
                    workbook = openpyxl.Workbook(write_only=False)
                    # Get the parent directories of each file
                    parent_directories = os.path.dirname(f)
                    # Remove characters before and including the last part of each directory
                    normalized_path = os.path.normpath(parent_directories)
                    modified_dir = os.path.basename(normalized_path)
                    counter = 1
                    while True:
                        folder = os.path.join(self.Save, f"{modified_dir}_Glyc_{counter}").replace("/", "\\")
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                            excel_file_path = f'{folder}\\output_Group{indx}_{ifl+1}_{modified_dir}.xlsx'
                            self.fo.append(folder)
                            sheet_name = modified_dir
                            # 在工作簿中创建一个新的工作表
                            sheet = workbook.create_sheet(sheet_name)
                            Peptide_Analysis = workbook.create_sheet('Peptide_Analysis')
                            # 写入表头
                            sheet.append(['Protein_ID', 'Peptide', 'Intensity'])
                            Intensity_Sum = 0
                            filtered_data = form[(form['total_signal'] != 0)]
                            glycopeptide = filtered_data['glycopeptide'].reset_index(drop=True)

                            total_signal = filtered_data['total_signal'].reset_index(drop=True)
                            protein_name = filtered_data['protein_name'].reset_index(drop=True)
                            total_signal_li = total_signal.tolist()
                            sites = filtered_data['n_glycosylation_sites'].reset_index(drop=True)
                            # 现为单一蛋白，2024.3.27，增添Protein ID当作识别方式，先循环获取，字典赋值
                            # 后面逐个输出
                            a = re.compile('\|[A-Za-z0-9]*\|')
                            protein_ID = []
                            for item in protein_name.tolist():
                                protein_ID.append(a.findall(item))
                            sites_list = sites.tolist()
                            glycopeptide_list = glycopeptide.tolist()
                            pattern = re.compile("\{.*\}")
                            # Glycan Ploting Data
                            for i,peptide in enumerate(glycopeptide_list):
                                glucose[i] = pattern.findall(peptide)
                                glucose_use_signal[i] = total_signal_li[i]
                            for item in total_signal_li:
                                Intensity_Sum += float(item)
                                # 逐行写入数据
                            for row_data in zip(filtered_data['protein_name'], form['glycopeptide'],form['total_signal']):
                                sheet.append(row_data)
                            R_Glycopeptide,NR_Glycopeptide = self.find_duplicates(glycopeptide_list)
                            tem_x = []
                            for ID in protein_ID:
                                tem_x.append(ID[0].split('|')[1])
                            protein_ID = tem_x
                            del tem_x
                            R_ID,NR_ID = self.find_duplicates(protein_ID)
                            # 此处要整改
                            for key,value in R_ID.items():
                                tem_x2 = []
                                tem_x1 = []
                                tem_y = []
                                tem_z = []
                                for value1 in value:
                                    tem_x2.append(sites_list[value1])
                                    tem_x1.append(total_signal_li[value1])
                                    tem_y.append(glucose[value1])
                                    tem_z.append(glucose_use_signal[value1])
                                sites_dic[ifl][key] = tem_x2
                                sites_signal_dic[ifl][key] = tem_x1
                                glycan_dic[ifl][key] = tem_y
                                glycan_signal_dic[ifl][key] = tem_z
                            for key,value in NR_ID.items():
                                sites_dic[ifl][key] = sites_list[value]
                                sites_signal_dic[ifl][key] = total_signal_li[value]
                                glycan_dic[ifl][key] = glucose[value]
                                glycan_signal_dic[ifl][key] = glucose_use_signal[value]
                            try:
                                del tem_x2, tem_y, tem_z, tem_x1
                            except Exception as e: # 检修用变量e
                                print('Without any repeat Proteins in process of one file')
                                pass
                            # 取对应的位点字典后，求强度，糖型相同操作，断点检查结合逐步运行
                            print(R_ID)
                            analysis_signal.write(R_ID)
                            logging.info(R_ID)
                            print(NR_ID)
                            analysis_signal.write(NR_ID)
                            logging.info(NR_ID)
                            default_sheet = workbook['Sheet']
                            workbook.remove(default_sheet)

                            Total_peptide = []
                            # 核心问题在于此
                            Rpart = {}
                            # Calculate the sum of all intensities for both repeating and non-repeating peptides
                            total_intensity_all = Intensity_Sum
                            self.both_intensity_all[indx][ifl] = total_intensity_all # 可用
                            # Loop through peptide_dealing_R
                            for key, values in R_Glycopeptide.items():
                                for value in values:
                                    ID = str(protein_ID[value])
                                    intensity = float(total_signal_li[value])
                                    # Use a tuple (ID, key) as a key for the dictionary
                                    key_tuple = (ID, key)
                                    # If the key is not in the dictionary, add it with the current intensity
                                    if key_tuple not in Rpart:
                                        Rpart[key_tuple] = intensity
                                    else:
                                        # If the key is already in the dictionary, update the intensity
                                        Rpart[key_tuple] += intensity
                            # Convert the dictionary to a list of lists
                            part_list1 = [[ID, key, intensity] for (ID, key), intensity in Rpart.items()]
                            Total_peptide.append(part_list1)
                            # print(part_list1)
                            # Write data to 'Peptide_Analysis' sheet for repeating peptides
                            Peptide_Analysis.append(['Protein Name', 'Glycopeptide', 'Total Intensity', 'Ratio'])
                            for row_data in part_list1:
                                # Calculate the ratio for each row
                                ratio = row_data[2] / total_intensity_all
                                Peptide_Analysis.append(row_data + [ratio])
                                # Loop through peptide_dealing_NR
                            part_list2 = []
                            for key, value in NR_Glycopeptide.items(): # 此处出问题
                                ID = str(protein_ID[value])
                                intensity = float(total_signal_li[value])
                                # Calculate the ratio for each non-repeating peptide
                                ratio = intensity / total_intensity_all
                                # Append data to 'part_list2' list for non-repeating peptides
                                part_list2.append([ID, key, intensity, ratio])
                                # Write data to 'Peptide_Analysis' sheet for non-repeating peptides
                            for row_data in part_list2:
                                Peptide_Analysis.append(row_data)
                            # print(part_list2)
                            Total_peptide.append(part_list2)
                            self.Both_Peptide[indx][ifl] = Total_peptide # 可用
                            analysis_signal.write(self.Both_Peptide)
                            logging.info(self.Both_Peptide)
                            # print(Both_Peptide)

                            # 保存Excel文件
                            # 设置自适应列宽(Sheet)
                            for column in sheet.columns:
                                max_length = 0
                                column = [cell for cell in column]
                                for cell in column:
                                    try:
                                        if len(str(cell.value)) > max_length:
                                            max_length = len(cell.value)
                                    except:
                                        pass
                                adjusted_width = (max_length + 2)
                                sheet.column_dimensions[column[0].column_letter].width = adjusted_width
                            # 设置自适应列宽(Peptide_Analysis)
                            for column in Peptide_Analysis.columns:
                                max_length = 0
                                column = [cell for cell in column]
                                for cell in column:
                                    try:
                                        if len(str(cell.value)) > max_length:
                                            max_length = len(cell.value)
                                    except:
                                        pass
                                adjusted_width = (max_length + 2)
                                Peptide_Analysis.column_dimensions[column[0].column_letter].width = adjusted_width
                            # print(folder)

                            workbook.save(excel_file_path)
                            workbook.close()
                            print(f"Excel File has saved at: {excel_file_path}")
                            analysis_signal.write(f"Excel File has saved at: {excel_file_path}")
                            logging.info(f"Excel File has saved at: {excel_file_path}")
                            break
                        counter += 1

                # 创建一个字典用于存储相同 n_glycosylation_sites 对应的 total_signal 总和，用于作图
                # Test Finished 20240421 23:07
                for key,value0 in sites_dic.items():
                    for key1,value1 in value0.items():
                        sum_by_nglyco = {}
                        value = value1
                        t = sites_signal_dic[key][key1]
                        if type(t) == float or type(t) == int:
                            total_signal_plot = t
                            if not re.findall(';', str(value)):
                                n_glycosylation_sites = str(value)
                                # 利用字典累加 total_signal
                                if n_glycosylation_sites in sum_by_nglyco:
                                    sum_by_nglyco[n_glycosylation_sites] += total_signal_plot
                                else:
                                    sum_by_nglyco[n_glycosylation_sites] = total_signal_plot
                            else:
                                l_site_ca = value.split(';')
                                for tem_site in l_site_ca:
                                    n_glycosylation_sites = str(tem_site)
                                    # 利用字典累加 total_signal
                                    if n_glycosylation_sites in sum_by_nglyco:
                                        sum_by_nglyco[n_glycosylation_sites] += total_signal_plot
                                    else:
                                        sum_by_nglyco[n_glycosylation_sites] = total_signal_plot
                        elif type(t) == list:
                            for g,total_signal_plot in enumerate(t):
                                site_ca = value[g]
                                if not re.findall(';', str(site_ca)):  # TypeError('expected string or bytes-like object')
                                    n_glycosylation_sites = str(site_ca)
                                    # 利用字典累加 total_signal
                                    if n_glycosylation_sites in sum_by_nglyco:
                                        sum_by_nglyco[n_glycosylation_sites] += total_signal_plot
                                    else:
                                        sum_by_nglyco[n_glycosylation_sites] = total_signal_plot
                                else:
                                    l_site_ca = site_ca.split(';')
                                    for tem_site in l_site_ca:
                                        n_glycosylation_sites = str(tem_site)
                                        # 利用字典累加 total_signal
                                        if n_glycosylation_sites in sum_by_nglyco:
                                            sum_by_nglyco[n_glycosylation_sites] += total_signal_plot
                                        else:
                                            sum_by_nglyco[n_glycosylation_sites] = total_signal_plot
                        new_sum_by_nglyco = {}
                        for inx,l in sum_by_nglyco.items():
                            if l != 0:
                                inx_new = int(inx)+1
                                l_new = l
                                new_sum_by_nglyco[inx_new] = l_new
                        new_sum_by_nglyco=sorted(new_sum_by_nglyco.items(), reverse=False)
                        del sum_by_nglyco
                        # 提取字典的键和值
                        nglyco_sites = []
                        total_signal_sum = []
                        for l in new_sum_by_nglyco:
                            nglyco_sites.append(str(l[0]))
                            total_signal_sum.append(l[1])

                        mpl.rcParams['font.family'] = 'Times New Roman'
                        # 创建一个新的图表
                        plt.figure(figsize=(10, 6))
                        # 绘制条形图
                        plt.bar(nglyco_sites, total_signal_sum,color='blue', alpha=0.7)
                        df = pd.DataFrame(({'Sites':nglyco_sites, 'Intens.':total_signal_sum}))

                        # x_maxsize = 200
                        # # change x internal size
                        # plt.gca().margins(x=0)
                        # plt.gcf().canvas.draw()
                        #
                        # # set size
                        # maxsize = x_maxsize
                        # m = 0.2
                        # N = len(nglyco_sites)
                        # s = maxsize / plt.gcf().dpi * N + 2 * m
                        # margin = m / plt.gcf().get_size_inches()[0]
                        #
                        # plt.gcf().subplots_adjust(left=margin, right=1. - margin)
                        # plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])

                        # 设置图表标题和坐标轴标签
                        plt.title(f'Total Signal Sum by N-Glycosylation Sites-Group{indx}-File{key+1}_{key1}')
                        plt.xlabel('N-Glycosylation Sites')
                        plt.ylabel('Total Signal Sum')
                        plt.autoscale(enable=True,axis='both',tight=None)
                        # Save
                        Save_Sites_Group = os.path.join(self.Save, f"Figure_GlycanSites_Group{indx}").replace("/", "\\")
                        if not os.path.exists(Save_Sites_Group):
                            os.makedirs(Save_Sites_Group)
                        plt.savefig(f"{Save_Sites_Group}\\{indx}_{key+1}_Glycan_sites_{key1}.png", dpi=1000)
                        plt.close()
                        df.to_excel(f"{Save_Sites_Group}\\{indx}_{key+1}_Glycan_sites_{key1}.xlsx", index=False)
                        analysis_signal.write(f"Glycosylation Site of {indx}_{key+1}_Glycan_sites_{key1} Figure Done!")

                # GlycanMass Plotting Module
                Hex = 180.16
                HexNAc = 221.21
                Fuc = 164.16
                Neu5Ac = 309.26
                Neu5Gc = 325.24
                for key0, value0 in glycan_dic.items():
                    glucose_plot_x = {}
                    glucose_plot_y_m = {}
                    x = {}
                    y = {}
                    for key1,value1 in value0.items():
                        glucose_plot_x[key1] = []
                        if type(value1) != list:
                            # 使用正则表达式提取键值对
                            item_1 = value1
                            matches = re.findall(r'(\w+|Neu5Ac|Neu5Gc):(\d+)', item_1)
                            # 创建一个列表存储单糖及其数量
                            sugar_counts = []
                            # 将提取的键值对添加到字典中
                            for match in matches:
                                sugar_name, count = match
                                sugar_counts.append((sugar_name,int(count)))
                            glucose_plot_x[key1].append(sugar_counts)
                        elif type(value1) == list:
                            glucose_plot_y_m[key1] = []
                            for i,item_1 in enumerate(value1):
                                if len(item_1) > 1:
                                    matches = re.findall(r'(\w+|Neu5Ac|Neu5Gc):(\d+)', item_1)
                                else:
                                    matches = re.findall(r'(\w+|Neu5Ac|Neu5Gc):(\d+)', item_1[0])
                                # 创建一个列表存储单糖及其数量
                                sugar_counts = []
                                # 将提取的键值对添加到字典中
                                for match in matches:
                                    sugar_name, count = match
                                    sugar_counts.append((sugar_name, int(count)))
                                glucose_plot_x[key1].append(sugar_counts)
                                glucose_plot_y_m[key1].append(i)
                    new_glucose = {}

                    for key, value in glucose_plot_x.items():
                        new_glucose[key] = []
                        for item_tem in value:
                            mass = 0
                            all_glycan_number= 0
                            for item_tem_tem in item_tem:
                                if item_tem_tem[0] == "Hex":
                                    mass += float(item_tem_tem[1]) * Hex
                                    all_glycan_number += item_tem_tem[1]
                                if item_tem_tem[0] == "HexNAc":
                                    mass += float(item_tem_tem[1]) * HexNAc
                                    all_glycan_number += item_tem_tem[1]
                                if item_tem_tem[0] == "Fuc":
                                    mass += float(item_tem_tem[1]) * Fuc
                                    all_glycan_number += item_tem_tem[1]
                                if item_tem_tem[0] == "Neu5Ac":
                                    mass += float(item_tem_tem[1]) * Neu5Ac
                                    all_glycan_number += item_tem_tem[1]
                                if item_tem_tem[0] == "Neu5Gc":
                                    mass += float(item_tem_tem[1]) * Neu5Gc
                                    all_glycan_number += item_tem_tem[1]
                            mass -= (all_glycan_number-1) * 18
                            mass -= 1
                            new_glucose[key].append(mass)

                    for key6,value6 in new_glucose.items():
                        x[key6] = value6
                    for key2, value2 in glucose_plot_y_m.items():
                        y[key2] = []
                        for item in value2:
                            t = glycan_signal_dic[key0][key2]
                            if type(t) == list:
                                y[key2].append(t[item])
                            else:
                                y[key2] = t
                    x_final = []
                    y_final = []
                    for key7, value7 in x.items():
                        for key8,value8 in enumerate(value7):
                            x_final.append(value8)
                            if type(y[key7]) == list:
                                y_final.append(y[key7][key8])
                            else:
                                y_final.append(y[key7])
                    del x,y
                    print(x_final)
                    analysis_signal.write(x_final)
                    print(y_final)
                    analysis_signal.write(y_final)
                    logging.info(f"Glycan-X Axis= {x_final}")
                    logging.info(f"Glycan-Y Axis= {y_final}")
                    x_mR,x_mNR = self.find_duplicates(x_final)
                    x_final_Calibration = []
                    y_final_Calibration = []
                    for x,index_list in x_mR.items():
                        x_final_Calibration.append(x)
                        y = 0
                        for inx in index_list:
                            y += y_final[inx]
                        y_final_Calibration.append(y)
                    for x,inx in x_mNR.items():
                        x_final_Calibration.append(x)
                        y_final_Calibration.append(y_final[inx])
                    print(x_final_Calibration)
                    print(y_final_Calibration)
                    analysis_signal.write(f"Glycan-X Axis= {x_final_Calibration}")
                    analysis_signal.write(f"Glycan-Y Axis= {y_final_Calibration}")
                    logging.info(f"Glycan-X Axis= {x_final_Calibration}")
                    logging.info(f"Glycan-Y Axis= {y_final_Calibration}")
                    mpl.rcParams['font.family'] = 'Times New Roman'
                    plt.rcParams['font.size'] = 12.5
                    # 创建一个新的图表
                    plt.figure(figsize=(10, 6))
                    # 绘制条形图
                    plt.bar(x_final_Calibration, y_final_Calibration, color='blue', alpha=0.7, width=3)
                    df = pd.DataFrame(({'Glycan Mass': x_final_Calibration, 'Intens.': y_final_Calibration}))
                    # x_maxsize = 200
                    # # change x internal size
                    # plt.gca().margins(x=0)
                    # plt.gcf().canvas.draw()
                    #
                    # # set size
                    # maxsize = x_maxsize
                    # m = 0.2
                    # N = len(x)
                    # s = maxsize / plt.gcf().dpi * N + 2 * m
                    # margin = m / plt.gcf().get_size_inches()[0]
                    #
                    # plt.gcf().subplots_adjust(left=margin, right=1. - margin)
                    # plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])

                    # 设置图表标题和坐标轴标签
                    plt.title(f'Total Signal of Glycosylation-Group{indx}-File{key0+1}')
                    plt.xlabel('N-Glycan Mass')
                    plt.ylabel('Total Signal Sum')
                    plt.autoscale(enable=True,axis='both',tight=None)
                    Save_Glycan_Group = os.path.join(self.Save, f"Figure_Glycan_Group{indx}").replace("/", "\\")
                    if not os.path.exists(Save_Glycan_Group):
                        os.makedirs(Save_Glycan_Group)
                    plt.savefig(f"{Save_Glycan_Group}\\{indx}_{key0+1}_Glycan.png", dpi=1000)
                    plt.close()
                    df.to_excel(f"{Save_Glycan_Group}\\{indx}_{key0+1}_Glycan.xlsx",index=False)
                    analysis_signal.write(f"{indx}_{key0+1} Glycan Total Signal Figure Done!")
                    logging.info(f"{indx}_{key0+1} Glycan Total Signal Figure Done!")
                glycan_final_dic = [glycan_dic,glycan_signal_dic]
                self.final[indx] = glycan_final_dic

            # the same group, different files
            for key0,value0 in self.Both_Peptide.items(): # Both_Peptide[index][ifl] = Total_peptide
                for key1,value1 in value0.items(): # value1 = ifl :  Total_peptide
                    ID_list = []
                    Peptide_list = []
                    Intensity_list = []
                    for value2 in value1:  # value2 = R_Peptide/NR_Peptide
                        for item in value2: # every Item in R_Peptide/NR_Peptide
                            tem_ID = item[0]
                            ID_list.append(tem_ID)
                            tem_Peptide = item[1]
                            Peptide_list.append(tem_Peptide)
                            tem_Intensity = item[2]
                            Intensity_list.append(tem_Intensity)
                    ID_list = tuple(ID_list)
                    Peptide_list = tuple(Peptide_list)
                    Intensity_list = tuple(Intensity_list)
                    # 输出重复元素与索引的函数部分，此处对元组进行操作
                    indices = {}
                    for index, item in enumerate(ID_list):
                        if item in indices:
                            indices[item].append(index)
                        else:
                            indices[item] = [index]
                    # print(ID_list)
                    # print(indices)

                    mpl.rcParams['font.family'] = 'Times New Roman'
                    x1_ID = []
                    y1_Intensity_Ratio = []
                    x2_Peptide = []
                    print(indices)
                    analysis_signal.write(indices)
                    logging.info(f'indices: {indices}')
                    for key in indices:
                        Intensity_Sum = 0
                        every_values = indices[key]
                        for item in every_values:
                            Intensity_Sum += Intensity_list[item]
                            x2_Peptide.append(Peptide_list[item])
                        a = re.search(r'[A-Za-z0-9]*', key).group(0) # IndexError('no such group')
                        x1_ID.append(a)
                        output = Intensity_Sum / self.both_intensity_all[key0][key1] # Glyc Fold-change运行关键
                        y1_Intensity_Ratio.append(output)
                    plt.rcParams['font.size'] = 12.5
                    fig, ax = plt.subplots(figsize=(18, 11))
                    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
                    bars = ax.bar(x1_ID, y1_Intensity_Ratio)
                    df = pd.DataFrame(({'Protein': x1_ID, 'Intens.': y1_Intensity_Ratio}))
                    self.Ratio3[(key0,key1)] = (x1_ID,y1_Intensity_Ratio)
                    # x_maxsize = 200
                    # # change x internal size
                    # plt.gca().margins(x=0)
                    # plt.gcf().canvas.draw()
                    #
                    # # set size
                    # maxsize = x_maxsize
                    # m = 0.2
                    # N = len(x1_ID)
                    # s = maxsize / plt.gcf().dpi * N + 2 * m
                    # margin = m / plt.gcf().get_size_inches()[0]
                    #
                    # plt.gcf().subplots_adjust(left=margin, right=1. - margin)
                    # plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])

                    for tick in ax.get_xticklabels():
                        tick.set_rotation(45)

                    # 使用循环遍历数据点
                    threshold = 0.1
                    for bar, value in zip(bars, y1_Intensity_Ratio):
                        # 根据条件判断是否添加注释
                        if value > threshold:
                            ax.annotate(f'{value * 100:.2f}%',  # 格式化数值
                                        xy=(bar.get_x() + bar.get_width() / 2, value),  # 获取柱形中心位置
                                        xytext=(0, 3),  # 调整文本位置
                                        textcoords="offset points",
                                        ha='center', va='bottom')
                    plt.title("Protein Ratio")
                    plt.autoscale(enable=True,axis='both',tight=None)
                    # plt.show()
                    Save_ProteinRatio_Group = os.path.join(self.Save, f"Figure_ProteinRatio_Group{key0}").replace("/", "\\")
                    if not os.path.exists(Save_ProteinRatio_Group):
                        os.makedirs(Save_ProteinRatio_Group)
                    plt.savefig(f"{Save_ProteinRatio_Group}\\{key0}_{key1 + 1}_ProteinRatio.png", dpi=1000)
                    plt.close()
                    df.to_excel(f"{Save_ProteinRatio_Group}\\{key0}_{key1 + 1}_ProteinRatio.xlsx", index=False)
                    analysis_signal.write(f"{key0}_{key1 + 1} Ratio Figure Done!")
                    logging.info(f"{key0}_{key1 + 1} Ratio Figure Done!")

            # the last part
            Ratio2 = []
            # print(Both_Peptide)
            tem1_3 = None
            tem2_3 = None
            print(len(self.Both_Peptide))
            analysis_signal.write(f"Length of Both_Prptide: {len(self.Both_Peptide)}")
            print(self.both_intensity_all)
            analysis_signal.write(f"Both_Prptide: {self.Both_Peptide}")
            for key0, value0 in self.Both_Peptide.items():  # 两组
                tem1_1 = {}  # 文件内的ID和
                tem2_1 = {}  # 文件内的ID和
                tem1_2 = {}  # 文件内的ID比值
                tem2_2 = {}  # 文件内的ID比值
                for key1, value1 in value0.items():  # 组内各个文件
                    ID_list = []
                    Peptide_list = []
                    Intensity_list = []
                    for value2 in value1:
                        for tem_list in value2:
                            tem_ID = tem_list[0]
                            ID_list.append(tem_ID)
                            tem_Peptide = tem_list[1]
                            Peptide_list.append(tem_Peptide)
                            tem_Intensity = tem_list[2]
                            Intensity_list.append(tem_Intensity)
                    R_ID, NR_ID = self.find_duplicates(ID_list)
                    for ID, item in R_ID.items():
                        ID_sum = 0
                        for number in item:
                            ID_sum += float(Intensity_list[number])
                        if key0 == 1:
                            tem1_1[(key1, ID)] = (ID, ID_sum)
                        else:
                            tem2_1[(key1, ID)] = (ID, ID_sum)
                    for ID, item in NR_ID.items():
                        ID_sum = 0
                        ID_sum += float(Intensity_list[item])
                        if key0 == 1:
                            tem1_1[(key1, ID)] = (ID, ID_sum)  # 问题在此
                        else:
                            tem2_1[(key1, ID)] = (ID, ID_sum)
                    for key123, value123 in self.both_intensity_all.items():  # 两组总峰强
                        if key0 == 1:
                            for key23, value23 in value123.items():  # 单个文件总峰强
                                for key124, value124 in tem1_1.items():  # 单个文件单一ID峰强
                                    if key124[0] == key23:  # key124为文件号，key23为单个文件总峰强字典文件号
                                        average = value124[1] / value23
                                        tem1_2[(key23, key124[1])] = (key124[1], average)
                        else:
                            for key23, value23 in value123.items():  # 单个文件总峰强
                                for key124, value124 in tem2_1.items():  # 单个文件单一ID峰强
                                    if key124[0] == key23:  # key124为文件号，key23为单个文件总峰强字典文件号
                                        average = value124[1] / value23
                                        tem2_2[(key23, key124[1])] = (key124[1], average)
                # 对3个文件的相同ID进行取平均，第一组
                # 用于存储相同第一个值的元组的第二个值的列表 20240625 0：22尚未测试完成
                value_dict = defaultdict(list)
                if key0 == 1:
                    for inner_key, (first, second) in tem1_2.items():
                        # 将第二个值添加到相同第一个值的列表中
                        value_dict[first].append(second)
                    # 计算平均值并保存到新字典中
                    average_dict = {key: sum(values) / len(values) for key, values in value_dict.items()}
                    tem1_3 = average_dict  # 单组内的ID平均值
                    del tem1_1, tem1_2
                # 对3个文件的相同ID进行取平均，第二组
                else:
                    for inner_key, (first, second) in tem2_2.items():
                        # 将第二个值添加到相同第一个值的列表中
                        value_dict[first].append(second)
                    # 计算平均值并保存到新字典中
                    average_dict = {key: sum(values) / len(values) for key, values in value_dict.items()}
                    tem2_3 = average_dict  # 单组内的ID平均值
                    del tem2_1, tem2_2
            for key1, value1 in tem1_3.items():
                for key2, value2 in tem2_3.items():
                    if value1 != 0 and value2 != 0 and key1 == key2:
                        output = value1 / value2
                        if log2(output) >= log2(self.Multiple_selection):  # 上调
                            Ratio2.append((key1, output))
                            self.Both_ID_list.append(key1)
                        if log2(output) <= log2(self.Multiple_selection):  # 下调 -0.585
                            Ratio2.append((key1, output))
                            self.Both_ID_list.append(key1)
            print(Ratio2)
            analysis_signal.write(Ratio2)
            logging.info(f'Ratio2: {Ratio2}')

            x1_ID = []
            y1_Intensity_Ratio = []
            # 2024.2.16 0:25 Ratio2_1.values(): > Ratio2_1
            for value1 in Ratio2:
                value = abs(value1[1])  # TypeError("bad operand type for abs(): 'str'")
                ID = re.search(r'[A-Z0-9]+', value1[0]).group(0)
                x1_ID.append(ID)
                y1_Intensity_Ratio.append(value)
            jug = len(Ratio2)
            i = 0
            # jugement to make a decision to deal with the data
            if jug != 0:
                for item in Ratio2:
                    if item[1] == 0:
                        i += 1
                if i != jug:
                    plt.rcParams['font.size'] = 12.5

                    fig, ax = plt.subplots(figsize=(18, 11))
                    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)

                    bars = ax.bar(x1_ID, y1_Intensity_Ratio)
                    df = pd.DataFrame(({'Protein': x1_ID, 'Intens.': y1_Intensity_Ratio}))
                    # x_maxsize = 200
                    # # change x internal size
                    # plt.gca().margins(x=0)
                    # plt.gcf().canvas.draw()
                    #
                    # # set size
                    # maxsize = x_maxsize
                    # m = 0.2
                    # N = len(x1_ID)
                    # s = maxsize / plt.gcf().dpi * N + 2 * m
                    # margin = m / plt.gcf().get_size_inches()[0]
                    #
                    # plt.gcf().subplots_adjust(left=margin, right=1. - margin)
                    # plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])

                    # ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=1))
                    # X轴角度45°
                    for tick in ax.get_xticklabels():
                        tick.set_rotation(45)

                    # 使用循环遍历数据点
                    threshold = 0.1
                    for bar, value in zip(bars, y1_Intensity_Ratio):
                        # 根据条件判断是否添加注释
                        if value > threshold:
                            ax.annotate('%.3f'% value,  # 格式化数值
                                        xy=(bar.get_x() + bar.get_width() / 2, value),  # 获取柱形中心位置
                                        xytext=(0, 3),  # 调整文本位置
                                        textcoords="offset points",
                                        ha='center', va='bottom')
                    plt.autoscale(True)
                    plt.title("Protein Ratio Difference Multiple")
                    plt.savefig(f"{self.Save}\\Protein_Ratio_Difference.png", dpi=1000)
                    plt.close()
                    df.to_excel(f"{self.Save}\\ProteinRatio.xlsx", index=False)
                    print("Done！")
                    analysis_signal.write("Finished Compared Part")
                    logging.info("Finished Compared Part")
                    # print("Please press Enter to exit~~")
                    # input()
                else:
                    print("Both sets of data are the same......")
                    analysis_signal.write("Both sets of data are the same......\nFinished Compared Part!")
                    logging.info("Both sets of data are the same......\nFinished Compared Part!")
                    # print("Please press Enter to exit~~")
                    # input()
            elif jug == 0:
                print("No relationship both of them......")
                analysis_signal.write("Both sets of data are the same\nFinished Compared Part!")
                logging.info("Both sets of data are the same\nFinished Compared Part!")
                # print("Please press Enter to exit~~")
                # input()
            try:
                list1 = set(self.Both_ID_list)
                list2 = list(list1)
                output_path = self.Save
                filename = 'NewProteinDatabase_After_Glyc'
                Total2 = list2
                fasta = self.fasta_input
                with open(file=os.path.join(output_path, filename + ".fasta"), mode="w", encoding="utf-8") as r:
                    for item1 in Total2:
                        for item2 in fasta:
                            op = re.search(re.escape(item1), item2)
                            if op:
                                matched_text = item2
                                output = matched_text
                                print(output)
                                r.write(output)
                    r.close()
            except Exception as er:
                logging.info(er)
                analysis_signal.write(er)
            print(self.Ratio3)
            return self.Ratio3,self.final
        except Exception as e:
            logging.info(e)
            analysis_signal.write(e)
        finally:
            print("Cleanup and exit.")

        # 恢复标准输出
        sys.stdout = sys.__stdout__
        analysis_signal.write("Compared Part of Analysis Thread completed.")
        logging.info("Compared Part of Analysis Thread completed.")