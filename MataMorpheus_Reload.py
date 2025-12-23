import threading

import pandas as pd
import csv
import os
import re
from AnalysisSignalModule import AnalysisSignal
outputSignal = AnalysisSignal()
class MetaMorpheus_Fasta_Reload(threading.Thread):
    @staticmethod
    def find_files(directory_path):
        tsv_files = None
        for root, dirs, files in os.walk(directory_path):
            if 'oglyco.psmtsv' in files:
                psmtsv_path = os.path.join(root, 'oglyco.psmtsv').replace("/", "\\")
                tsv_files = psmtsv_path
                print(f"Found the psmtsv file: {psmtsv_path}")
        #         outputSignal.write(psmtsv_path)
        # outputSignal.write(tsv_files)
        return tsv_files
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

    def __init__(self,psm1_path,psm2_path,fasta_file,save_path):
        super().__init__()
        self.total_output = []
        self.d = [psm1_path, psm2_path]
        self.both_ID = []
        self.fasta_file = fasta_file

    def run(self):
        for p,psm_path in enumerate(self.d):
            if type(psm_path) == str:
                path = self.find_files(psm_path)
                outputSignal.write(path)
                # self.update_text_browser1(f"{path}")
                # parent_directories = (os.path.dirname(path))
                # Remove characters before and including the last part of each directory
                # modified_directories = parent_directories.split("\\")[-1]
                # modified_dir = modified_directories.split("/")[-1]

                Total1 = []
                i = os.path.join(path)
                print(i)
                outputSignal.write(i)
                protein_ID = []
                with open(i, 'r') as file:
                    # 创建 CSV 读取器
                    reader = csv.reader(file, delimiter='\t')
                    header = next(reader)
                    # # 遍历每一行数据
                    for row in reader:
                        protein_ID.append(row[7])
                print(protein_ID)
                outputSignal.write(protein_ID)

                for match in protein_ID:
                    Total1.append(match)

                Total2 = set(Total1)
                Total2 = list(Total2)
                # print(Total2)
                for i in Total2:
                    self.both_ID.append(i)
                # outputSignal.write(Total2)
            elif type(psm_path) == list:
                for item in psm_path:
                    path = self.find_files(item)
                    outputSignal.write(path)
                    Total1 = []
                    i = os.path.join(path)
                    print(i)
                    outputSignal.write(i)
                    protein_ID = []
                    with open(i, 'r') as file:
                        # 创建 CSV 读取器
                        reader = csv.reader(file, delimiter='\t')
                        header = next(reader)
                        # # 遍历每一行数据
                        for row in reader:
                            protein_ID.append(row[7])
                    print(protein_ID)
                    outputSignal.write(protein_ID)

                    for match in protein_ID:
                        Total1.append(match)

                    Total2 = set(Total1)
                    Total2 = list(Total2)
                    # print(Total2)
                    for i in Total2:
                        self.both_ID.append(i)
                    # outputSignal.write(Total2)
        del Total2, Total1
        fasta_path = self.fasta_file
        tem_i = set(self.both_ID)
        Total2 = list(tem_i)
        with open(file=fasta_path, mode="r", encoding="utf-8") as o:
            tem = o.read()
            fasta = tem.split(sep=">")
            # with open(file=os.path.join(output_path,filename+".fasta"), mode="w", encoding="utf-8") as r:
            for item1 in Total2:
                for item2 in fasta:
                    a = re.compile(item1)
                    op = a.findall(item2)
                    for item3 in op:
                        if len(op) != 0:
                            sy = MetaMorpheus_Fasta_Reload.find_index_with_regex(fasta, item3)
                            output = ">" + fasta[sy]
                            print(output)
                            # r.write(output)
                            self.total_output.append(output)
                            # outputSignal.write(output)
            # r.close()
        o.close()
        return self.total_output
        # print("Done!")
        outputSignal.write("Done!")
        outputSignal.write('Fasta File Reload Done!')