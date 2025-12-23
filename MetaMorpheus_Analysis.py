# Mode: After MetaMorpheus
import os,logging,openpyxl,sys,csv
import re,pickle,pathlib
import pandas as pd
from AnalysisSignalModule import AnalysisSignal
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
from MetaMorpheus_FlashLFQ_Intensity import MetaMorpheus_FlashLFQ_Intensity
from math import log,log2
from collections import defaultdict
from resource_path import resource_path
MetaMorpheus_output = AnalysisSignal()
class Analysis_MetaMorpheus():
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
    def get_key_by_min_value(dictionary):
        fixed_index = 0
        purpose_index = 1
        min_val_k = []
        min_val = []
        for sub_list in dictionary.values():
            min_val.append(sub_list[fixed_index])
            min_val_k.append(sub_list[purpose_index])
        tem = min(min_val)
        for index,value in enumerate(min_val):
            if value == tem:
                return min_val_k[index]

    @staticmethod
    def qvalue_filter(qvalue,list):
        tem = []
        for j,element in enumerate(qvalue):
            if float(element) <= 0.05:
                tem.append(list[j])
        return tem

    @staticmethod
    def removeDuplicates(list):
        if len(list) == 0: return 0
        k = 1
        for i in range(1, len(list)):
            if list[i] != list[i - 1]:
                list[k] = list[i]
                k += 1
        return k
    @staticmethod
    def find_QuantifiedPeaks(directory_path):
        psmtsv_files = None
        for root, dirs, files in os.walk(directory_path):
            if 'QuantifiedPeaks.tsv' in files:
                psmtsv_path = os.path.join(root, 'QuantifiedPeaks.tsv').replace("/", "\\")
                psmtsv_files = psmtsv_path
                print(f"Found the tsv file: {psmtsv_path}")
        #         MetaMorpheus_output.write(psmtsv_path)
        # MetaMorpheus_output.write(psmtsv_files)
        return psmtsv_files

    @staticmethod
    def parse_index_range(index_range_str):
        # Extract the two numbers from the string '[117 to 126]' and return them as a tuple
        indices = re.findall(r'\d+', index_range_str)
        return int(indices[0]), int(indices[1])
    @staticmethod
    def find_o_glycosylation_sites(seq, index_range_strings):
        # index_range_strings = ''
        if seq != None and index_range_strings != '':
            o_glyco_sites = []
            # Parse the index range string
            start_index, end_index = Analysis_MetaMorpheus.parse_index_range(index_range_strings)
            # Remove bracketed annotations from the sequence (e.g., "[O-Glycosylation:H1N1G2 on X]")
            cleaned_seq = re.sub(r'\[.*?\]', '', seq)

            # Find positions of O-glycosylation in the original sequence
            o_glyco_pos = [(m.start(), m.group()[0]) for m in re.finditer(r'[ST]\[O-Glycosylation.*?\]', seq)]

            # Calculate the corrected positions after stripping the annotations
            corrected_positions = []
            for pos, aa in o_glyco_pos:
                # Compute the number of characters skipped (annotations) before this position
                annotations_before_pos = re.findall(r'\[.*?\]', seq[:pos])
                shift = sum(len(ann) for ann in annotations_before_pos)  # Sum of bracketed parts' lengths

                # Adjust the position based on how much we've stripped and the provided index range
                corrected_pos = pos - shift + start_index
                corrected_positions.append((aa, corrected_pos))  # (amino acid, corrected position)

            o_glyco_sites.extend(corrected_positions)  # Use extend to add all items at once

            return o_glyco_sites
        else:
            return None

    # @staticmethod
    # def find_oglyco_psmtsv(directory_path):
    #     psmtsv_files = []
    #     for root, dirs, files in os.walk(directory_path):
    #         if 'oglyco.psmtsv' in files:
    #             psmtsv_path = os.path.join(root, 'oglyco.psmtsv').replace("/", "\\")
    #             psmtsv_files.append(psmtsv_path)
    #             print(f"Found the psmtsv file: {psmtsv_path}")
    #             MetaMorpheus_output.write(psmtsv_path)
    #     MetaMorpheus_output.write(psmtsv_files)
    #     return psmtsv_files
    # @staticmethod
    # def find_olocationsites_file(directory_path):
    #     tsv_files = []
    #     directory_path = directory_path[0]
    #     for root, dirs, files in os.walk(directory_path):
    #         if 'protein_oglyco_localization.tsv' in files:
    #             psmtsv_path = os.path.join(root, 'protein_oglyco_localization.tsv').replace("/", "\\")
    #             tsv_files.append(psmtsv_path)
    #             print(f"Found the tsv file: {psmtsv_path}")
    #             MetaMorpheus_output.write(psmtsv_path)
    #     MetaMorpheus_output.write(tsv_files)
    #     return tsv_files

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
    def __init__(self, file_psm1, file_psm2, peak1, peak2, Save_path, fasta_input, Multiple_selection):

        self.path = {1:file_psm1, 2:file_psm2}
        self.Peak_path = {1:peak1,2:peak2}
        self.Save = Save_path
        self.fo = []
        self.Both_Peptide = {}
        self.both_intensity_all = {}
        self.Both_ID_list = []
        self.fasta_input = fasta_input
        self.Multiple_selection = Multiple_selection
        self.Ratio3 = {}
        MetaMorpheus_output.write(text=f"Starting analysis\nPath Group1: {file_psm1}, \nPath Group2: {file_psm2}, \nSave path: {Save_path},\nRAW Group1: {peak1},\nRAW Group2: {peak2}")
        logging.info(f"Starting analysis with \nPath1: {file_psm1}, \nPath2: {file_psm2}, \nSave path: {Save_path},\nRAW1: {peak1},\nRAW2: {peak2}")

    def run(self):
        try:
            self.final = {}
            for INDX0,item in self.path.items():
                self.final[INDX0] = {}
                glycan_dic = {}
                glycan_intensity_dic = {}
                self.Both_Peptide[INDX0] = {}
                self.both_intensity_all[INDX0] = {}
                glycan_for_final = {}
                glycan_slgnal_for_final = {}
                for idl,item0 in enumerate(item):
                    glycan_for_final[idl] = {}
                    glycan_slgnal_for_final[idl] = {}
                    i = item0  # PSM1
                    j = self.Peak_path[INDX0][idl]  # PEAKS1
                    # i = self.find_oglyco_psmtsv(i)[0]
                    PSM_DataFrame = MetaMorpheus_FlashLFQ_Intensity.read_sequences_from_csv(i)
                    Peak_DataFrame = MetaMorpheus_FlashLFQ_Intensity.read_sequences_from_tsv(j)
                    try:
                        container = MetaMorpheus_FlashLFQ_Intensity.run(PSM_DataFrame, Peak_DataFrame)
                        exit_code = container.wait()
                        logging.info(f"Container {container.id} finished with exit code {exit_code}")
                        if exit_code['StatusCode'] == 0:
                            logging.info("Command executed successfully, proceeding with next steps or files...")
                            print('Command executed successfully, proceeding with next steps or files...')
                        else:
                            logging.error(f"Command execution failed, check the container logs for details. Return exit_code: {exit_code['StatusCode']}")
                    except Exception as e:
                        print(e)
                    try:
                        with open(pathlib.Path(resource_path('resources')) / 'task_Meta' / 'cache_file_task_Meta.pkl', 'rb') as f:
                            DataDf = pickle.load(f)
                    except FileNotFoundError as e:
                        print(e)
                    # 'ID': peak_row['Protein Group'],
                    # 'Peptide': peak_row['Full Sequence'],
                    # 'Glycan': psm_row['Plausible GlycanComposition'],
                    # 'GlycanMass': psm_row['GlycanMass'],
                    # 'Intensity': peak_row['Peak intensity'],
                    # 'Site': find_o_glycosylation_sites(peak_row['Full Sequence'],
                    #                                    psm_row['Start and End Residues In Protein'])

                    # 表头如上，24.12.28 21：14
                    b = re.compile(r"DECOY_[A-Za-z0-9]*")
                    Protein_ID = DataDf['ID'].tolist()
                    Peptide = DataDf['Peptide'].tolist()
                    GlycanMass = DataDf['GlycanMass'].tolist()
                    GlycanStruc = DataDf['Glycan'].tolist()
                    intensity_after = DataDf['Intensity'].tolist()
                    Site = DataDf['Site'].tolist()
                    print(intensity_after)
                    logging.info(intensity_after)
                    MetaMorpheus_output.write(intensity_after)
                    R_ID, NR_ID = self.find_duplicates(Protein_ID)
                    print("Repeat Protein:",R_ID)
                    sites_peptide_signal_dic = {}
                    glycanmass_peptide_dic = {}
                    glycan_peptide_signal_dic = {}
                    site_n = {}
                    # 20240510 22：11
                    for key, value in R_ID.items():
                        tem_x1 = []
                        tem_x2 = []
                        tem_x3 = []
                        tem_x4 = []
                        tem_x5 = []
                        for value1 in value:
                            tem_x1.append(intensity_after[value1])
                            tem_x2.append(GlycanStruc[value1])
                            tem_x3.append(intensity_after[value1])
                            tem_x4.append(GlycanMass[value1])
                            tem_x5.append(Site[value1])
                        sites_peptide_signal_dic[key] = tem_x1
                        glycan_dic[key] = tem_x2
                        glycan_peptide_signal_dic[key] = tem_x3

                        # For final_analysis module
                        glycan_for_final[idl][key] = tem_x2
                        glycan_slgnal_for_final[idl][key] = tem_x3

                        glycanmass_peptide_dic[key] = tem_x4
                        site_n[key] = tem_x5
                        del tem_x1, tem_x2, tem_x3, tem_x4, tem_x5
                    for key, value in NR_ID.items():
                        site_n[key] = [Site[value]]
                        sites_peptide_signal_dic[key] = [intensity_after[value]]
                        glycan_dic[key] = GlycanStruc[value]
                        glycan_peptide_signal_dic[key] = intensity_after[value]
                        glycanmass_peptide_dic[key] = GlycanMass[value]

                        # For final_analysis module
                        glycan_for_final[idl][key] = GlycanStruc[value]
                        glycan_slgnal_for_final[idl][key] = intensity_after[value]

                    self.final[INDX0] = [glycan_for_final, glycan_slgnal_for_final]  # Fix 250418
                    workbook = openpyxl.Workbook(write_only=False)
                    # Get the parent directories of each file
                    parent_directories = (os.path.dirname(item0)) # idl is an int type data
                    # Remove characters before and including the last part of each directory
                    modified_directories = parent_directories.split("/")[-1]
                    modified_dir = modified_directories

                    counter = 1
                    while True:
                        folder = os.path.join(self.Save, f"{modified_dir}_O-Pair_{counter}").replace("/", "\\")
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                            excel_file_path = f'{folder}\\output_Group{INDX0}_File{idl+1}_{modified_dir}.xlsx'
                            self.fo.append(folder)
                            sheet_name = modified_dir
                            # 在工作簿中创建一个新的工作表
                            sheet = workbook.create_sheet(sheet_name)
                            Peptide_Analysis = workbook.create_sheet('Peptide_Analysis')
                            # 写入表头
                            sheet.append(['Protein_ID', 'Peptide', 'Intensity'])
                            total_signal = intensity_after
                            # both_intensity_all.append(total_signal)
                            protein_name = Protein_ID
                            # 逐行写入数据
                            for row_data in zip(Protein_ID, Peptide, intensity_after):
                                sheet.append(row_data)
                            R_Glycopeptide, NR_Glycopeptide = self.find_duplicates(Peptide)
                            default_sheet = workbook['Sheet']
                            workbook.remove(default_sheet)
                            Total_peptide = []
                            # 核心问题在于此
                            Rpart = {}
                            # Calculate the sum of all intensities for both repeating and non-repeating peptides
                            total_intensity_all = sum(float(value) for value in total_signal)
                            self.both_intensity_all[INDX0][idl] = total_intensity_all
                            # Loop through peptide_dealing_R
                            for key, values in R_Glycopeptide.items():
                                for value in values:
                                    ID = str(protein_name[value])
                                    intensity = float(total_signal[value])
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
                            print(NR_Glycopeptide)
                            print(protein_name)
                            print(total_signal)
                            for key, value in NR_Glycopeptide.items():
                                ID = str(protein_name[value])
                                intensity = float(total_signal[value])
                                # Calculate the ratio for each non-repeating peptide
                                ratio = intensity / total_intensity_all
                                # Append data to 'part_list2' list for non-repeating peptides
                                part_list2.append([ID, key, intensity, ratio])
                                # Write data to 'Peptide_Analysis' sheet for non-repeating peptides
                            for row_data in part_list2:
                                Peptide_Analysis.append(row_data)
                            # print(part_list2)
                            Total_peptide.append(part_list2)
                            self.Both_Peptide[INDX0][idl] = Total_peptide
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
                            MetaMorpheus_output.write(text=f"Excel File has saved at: {excel_file_path}")
                            logging.info(f"Excel File has saved at: {excel_file_path}")
                            break
                        counter += 1

                    # 创建一个字典用于存储相同o_glycosylation_sites对应的 total_signal 总和，用于作图
                    # 着重分析清楚，site_ID_final的数据组成结构，来实现功能 25.1.1
                    try:
                        b = re.compile(r"DECOY_[A-Za-z0-9]*")
                        # Key: ID, value: Site
                        for key, value in site_n.items():
                            sum_by_oglyco = {}
                            if not b.findall(key) and isinstance(value, list):
                                for idex, value1 in enumerate(value):
                                    total_signal_plot = float(sites_peptide_signal_dic[key][idex])
                                    if isinstance(value1, list):
                                        for l1 in value1:
                                            o_glycosylation_sites = str(l1)
                                            sum_by_oglyco.setdefault(o_glycosylation_sites, 0)
                                            sum_by_oglyco[o_glycosylation_sites] += total_signal_plot
                                    else:
                                        o_glycosylation_sites = str(value1)
                                        sum_by_oglyco.setdefault(o_glycosylation_sites, 0)
                                        sum_by_oglyco[o_glycosylation_sites] += total_signal_plot
                            if sum_by_oglyco:
                                # 转换键为整数，排序，然后转换回字符串
                                # 创建一个空列表用于存储转换后的元组
                                sorted_sum_by_oglyco = []
                                # 遍历字典的键值对
                                for k, v in sum_by_oglyco.items():
                                    # 使用正则表达式从键中提取数字部分
                                    match = re.search(r"\d+", k)
                                    if match:
                                        # 将提取的数字转换为整数，并作为元组添加到列表中
                                        sorted_sum_by_oglyco.append((int(match.group()), v))
                                # 对列表中的元组按照第一个元素（提取的整数）进行排序
                                sorted_sum_by_oglyco = sorted(sorted_sum_by_oglyco)
                                oglyco_sites, total_signal_sum = zip(*sorted_sum_by_oglyco)
                                oglyco_sites = [str(site) for site in oglyco_sites]  # 转换回字符串用于标签

                                plt.figure(figsize=(10, 6))
                                ax = plt.gca()
                                ax.bar(range(len(oglyco_sites)), total_signal_sum, color='blue', alpha=0.7)
                                df = pd.DataFrame(({'Sites': oglyco_sites, 'Intens.': total_signal_sum}))
                                ax.set_title(f'Total Signal Sum by O-Glycosylation Sites-Group{INDX0}-File{idl + 1}-{key}')
                                ax.set_xlabel('O-Glycosylation Sites')
                                ax.set_ylabel('Total Signal Sum')

                                # 根据标签数量和图表宽度调整显示的标签数量
                                total_labels = len(oglyco_sites)
                                step = max(1, total_labels // 15)  # For example, a label is displayed every 15 points, which can be adjusted as needed
                                plt.xticks(range(len(oglyco_sites))[::step], oglyco_sites[::step])

                                plt.tight_layout()
                                save_path = os.path.join(self.Save, f"Figure_GlycanSites_Group{INDX0}")
                                os.makedirs(save_path, exist_ok=True)
                                k = re.sub(r'[^\w]', '-', key)
                                plt.savefig(os.path.join(save_path, f"{INDX0}_{idl + 1}_Glycan_sites_{k}.png"), dpi=1000)
                                x_p = os.path.join(save_path, f"{INDX0}_{idl + 1}_Glycan_sites_{k}.xlsx")
                                df.to_excel(x_p,index=False)
                                plt.close()
                    except Exception as e:
                        print(e)
                        MetaMorpheus_output.write(e)


                    # GlycanMass Plotting Module
                    try:
                        sum_by_oglyco_glycan = {}
                        for key,value in glycanmass_peptide_dic.items():
                            if b.findall(key) == [] and isinstance(value, list):
                                for idex,value1 in enumerate(value):
                                    total_signal_plot = float(glycan_peptide_signal_dic[key][idex])
                                    o_glycosylation_glycan = float(value1)
                                    # 利用字典累加 total_signal
                                    if o_glycosylation_glycan in sum_by_oglyco_glycan:
                                        sum_by_oglyco_glycan[o_glycosylation_glycan] += total_signal_plot
                                    else:
                                        sum_by_oglyco_glycan[o_glycosylation_glycan] = total_signal_plot
                            elif b.findall(key) == [] and isinstance(value, (int, float, str)):
                                total_signal_plot = float(glycan_peptide_signal_dic[key])
                                o_glycosylation_glycan = float(value)
                                # 利用字典累加 total_signal
                                if o_glycosylation_glycan in sum_by_oglyco_glycan:
                                    sum_by_oglyco_glycan[o_glycosylation_glycan] += total_signal_plot
                                else:
                                    sum_by_oglyco_glycan[o_glycosylation_glycan] = total_signal_plot
                            else:
                                pass
                        print(sum_by_oglyco_glycan)
                        x = []
                        y = []
                        # 提取字典的键和值
                        for keym,valuem in sum_by_oglyco_glycan.items():
                            x.append(float(keym))
                            y.append(float(valuem))
                        mpl.rcParams['font.family'] = 'Times New Roman'
                        plt.rcParams['font.size'] = 12.5
                        # 创建一个新的图表
                        plt.figure(figsize=(10, 6))
                        # 绘制条形图
                        plt.bar(x, y, color='blue', alpha=0.7, width=3)
                        df = pd.DataFrame(({'Glycan Mass': x, 'Intens.': y}))
                        # 设置图表标题和坐标轴标签
                        plt.title(f'Total Signal of Glycosylation-Group{INDX0}-File{idl+1}')
                        plt.xlabel('O-Glycan Mass')
                        plt.ylabel('Total Signal Sum')
                        plt.autoscale()
                        Save_Glycan_Group = os.path.join(self.Save, f"Figure_Glycan_Group{INDX0}").replace("/","\\")
                        if not os.path.exists(Save_Glycan_Group):
                            os.makedirs(Save_Glycan_Group)
                        plt.savefig(f"{Save_Glycan_Group}\\{INDX0}_{idl + 1}_Glycan.png", dpi=1000)
                        plt.close()
                        df.to_excel(f"{Save_Glycan_Group}\\{INDX0}_{idl + 1}_Glycan.xlsx",index=False)
                        MetaMorpheus_output.write("Glycan Total Signal Figure Done!")
                        logging.info(f"{INDX0}_{idl + 1} Glycan Total Signal Figure Done!")
                    except Exception as e:
                        print(e)
                        MetaMorpheus_output.write(e)
                    # glycan_final_dic = Remember the logic to complete the dict ^_^

            # the same group, different files
            for key0, value0 in self.Both_Peptide.items():  # Both_Peptide[index][ifl] = Total_peptide
                for key1,value1 in value0.items(): # value1 = ifl :  Total_peptide
                    ID_list = []
                    Peptide_list = []
                    Intensity_list = []
                    for value2 in value1: # value2 = R_Peptide/NR_Peptide
                        for tem_list in value2:   # every Item in R_Peptide/NR_Peptide
                            tem_ID = tem_list[0]
                            ID_list.append(tem_ID)
                            tem_Peptide = tem_list[1]
                            Peptide_list.append(tem_Peptide)
                            tem_Intensity = tem_list[2]
                            Intensity_list.append(tem_Intensity)
                    ID_list = tuple(ID_list)
                    Peptide_list = tuple(Peptide_list)
                    Intensity_list = tuple(Intensity_list)
                    # 输出重复元素与索引的函数部分，此处对元组进行操作
                    indices = {}
                    for indx, item in enumerate(ID_list):
                        if item in indices:
                            indices[item].append(indx)
                        else:
                            indices[item] = [indx]
                    # print(ID_list)
                    # print(indices)

                    mpl.rcParams['font.family'] = 'Times New Roman'
                    x1_ID = []
                    y1_Intensity_Ratio = []
                    x2_Peptide = []
                    for key in indices:
                        Intensity_Sum = 0
                        every_values = indices[key]
                        for item in every_values:
                            Intensity_Sum += Intensity_list[item]
                            x2_Peptide.append(Peptide_list[item])
                        x1_ID.append(key)
                        output = Intensity_Sum / self.both_intensity_all[key0][key1]
                        y1_Intensity_Ratio.append(output)
                    plt.rcParams['font.size'] = 12.5
                    fig, ax = plt.subplots(figsize=(18, 11))
                    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
                    bars = ax.bar(x1_ID, y1_Intensity_Ratio)
                    df = pd.DataFrame(({'Protein': x1_ID, 'Intens.': y1_Intensity_Ratio}))
                    self.Ratio3[(key0,key1)] = (x1_ID,y1_Intensity_Ratio)
                    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=1))
                    total_labels = len(x1_ID)
                    step = max(1,total_labels // 15)  # For example, a label is displayed every 15 points, which can be adjusted as needed
                    plt.xticks(range(len(x1_ID))[::step], x1_ID[::step])
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
                    plt.autoscale()
                    Save_ProteinRatio_Group = os.path.join(self.Save, f"Figure_ProteinRatio_Group{key0}").replace("/", "\\")
                    if not os.path.exists(Save_ProteinRatio_Group):
                        os.makedirs(Save_ProteinRatio_Group)
                    plt.savefig(f"{Save_ProteinRatio_Group}\\{key0}_{key1 + 1}_ProteinRatio.png", dpi=1000)
                    plt.close()
                    df.to_excel(f"{Save_ProteinRatio_Group}\\{key0}_{key1 + 1}_ProteinRatio.xlsx",index=False)
                    MetaMorpheus_output.write(f"{key0}_{key1 + 1} Ratio Figure Done!")
                    logging.info(f"{key0}_{key1 + 1} Ratio Figure Done!")

            Ratio2 = []
            # print(Both_Peptide)
            tem1_3 = None
            tem2_3 = None
            print(len(self.Both_Peptide))
            MetaMorpheus_output.write(f"Length of Both_Prptide: {len(self.Both_Peptide)}")
            print(self.both_intensity_all)
            MetaMorpheus_output.write(f"Both_Prptide: {self.Both_Peptide}")
            for key0, value0 in self.Both_Peptide.items(): # 两组
                tem1_1 = {} # 文件内的ID和
                tem2_1 = {} # 文件内的ID和
                tem1_2 = {} # 文件内的ID比值
                tem2_2 = {} # 文件内的ID比值
                for key1,value1 in value0.items(): # 组内各个文件
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
                    R_ID,NR_ID = self.find_duplicates(ID_list)
                    for ID,item in R_ID.items():
                        ID_sum = 0
                        for number in item:
                            ID_sum += float(Intensity_list[number])
                        if key0 == 1:
                            tem1_1[(key1,ID)] = (ID, ID_sum)
                        else:
                            tem2_1[(key1,ID)] = (ID, ID_sum)
                    for ID,item in NR_ID.items():
                        ID_sum = 0
                        ID_sum += float(Intensity_list[item])
                        if key0 == 1:
                            tem1_1[(key1, ID)] = (ID, ID_sum)  # 问题在此
                        else:
                            tem2_1[(key1, ID)] = (ID, ID_sum)

                    for key123,value123 in self.both_intensity_all.items(): #两组总峰强
                        if key0 == 1:
                            for key23,value23 in value123.items():  # 单个文件总峰强
                                for key124,value124 in tem1_1.items():  # 单个文件单一ID峰强
                                    if key124[0] == key23:  # key124为文件号，key23为单个文件总峰强字典文件号
                                        average = value124[1]/value23
                                        tem1_2[(key23,key124[1])] = (key124[1],average)
                        else:
                            for key23,value23 in value123.items():  # 单个文件总峰强
                                for key124,value124 in tem2_1.items():  # 单个文件单一ID峰强
                                    if key124[0] == key23:  # key124为文件号，key23为单个文件总峰强字典文件号
                                        average = value124[1]/value23
                                        tem2_2[(key23,key124[1])] = (key124[1],average)
                # 对3个文件的相同ID进行取平均，第一组
                # 用于存储相同第一个值的元组的第二个值的列表 20240625 0：22尚未测试完成
                value_dict = defaultdict(list)
                if key0 == 1:
                    for inner_key, (first, second) in tem1_2.items():
                        # 将第二个值添加到相同第一个值的列表中
                        value_dict[first].append(second)
                    # 计算平均值并保存到新字典中
                    average_dict = {key: sum(values) / len(values) for key, values in value_dict.items()}
                    tem1_3 = average_dict # 单组内的ID平均值
                    del tem1_1,tem1_2
                # 对3个文件的相同ID进行取平均，第二组
                else:
                    for inner_key, (first, second) in tem2_2.items():
                        # 将第二个值添加到相同第一个值的列表中
                        value_dict[first].append(second)
                    # 计算平均值并保存到新字典中
                    average_dict = {key: sum(values) / len(values) for key, values in value_dict.items()}
                    tem2_3 = average_dict # 单组内的ID平均值
                    del tem2_1, tem2_2
            for key1, value1 in tem1_3.items():
                for key2, value2 in tem2_3.items():
                    if value1 != 0 and value2 != 0 and key1 == key2:
                        output = value1 / value2
                        if log2(output) >= log2(self.Multiple_selection):  # 上调
                            Ratio2.append((key1, output))
                            self.Both_ID_list.append(key1)
                        if log2(output) <= -log2(self.Multiple_selection):  # 下调 -0.585
                            Ratio2.append((key1, output))
                            self.Both_ID_list.append(key1)
            print(Ratio2)
            # MetaMorpheus_output.write(Ratio2)
            logging.info(f"Ratio2: {Ratio2}")

            x1_ID = []
            y1_Intensity_Ratio = []
            # 2024.2.16 0:25 Ratio2_1.values(): > Ratio2_1
            for value1 in Ratio2:
                value = abs(value1[1])
                x1_ID.append(value1[0])
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
                    # ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=1))
                    # X轴角度45°
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
                    plt.title("Protein Ratio Different Value")
                    plt.autoscale()
                    plt.savefig(f"{self.Save}\\Protein_Ratio_Difference.png", dpi=1000)
                    plt.close()
                    df.to_excel(f"{self.Save}\\Protein_Ratio_Difference.xlsx", index=False)
                    print("Done！")
                    MetaMorpheus_output.write(text="Finished Compared Part")
                    logging.info("Finished Compared Part")
                    # print("Please press Enter to exit~~")
                    # input()
                else:
                    print("Both sets of data are the same......")
                    MetaMorpheus_output.write(text="Both sets of data are the same......\nFinished Compared Part!")
                    logging.info("Both sets of data are the same......\nFinished Compared Part!")
                    # print("Please press Enter to exit~~")
                    # input()
            elif jug == 0:
                print("No relationship both of them......")
                MetaMorpheus_output.write(text="Both sets of data are the same\nFinished Compared Part!")
                logging.info("Both sets of data are the same\nFinished Compared Part!")
                # print("Please press Enter to exit~~")
                # input()
            # ERROR:root:An error occurred during analysis: 'list' object has no attribute 'items' 20240509 22:17
            try:
                list1 = set(self.Both_ID_list)
                list2 = list(list1)
                output_path = self.Save
                filename = 'NewProteinDatabase_After_O-Pair'
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
                MetaMorpheus_output.write(er)
            return self.Ratio3,self.final
        except Exception as e:
            logging.info(e)
            MetaMorpheus_output.write(e)
        finally:
            print("Cleanup and exit.")
    # 恢复标准输出
    sys.stdout = sys.__stdout__
    MetaMorpheus_output.write(text="Compared Part of Analysis Thread completed.")
    logging.info("Compared Part of Analysis Thread completed.")