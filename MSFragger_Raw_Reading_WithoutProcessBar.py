# Mode: After MSFragger
from pymsfilereader import MSFileReader
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import os
from AnalysisSignalModule import AnalysisSignal

AnalysisSignal = AnalysisSignal()
import logging

# logging.level(message)创建一条level级别的日志
# logging.debug("This is a debug log")
# logging.info("This is a info log")
# logging.warning("This is a warning log")
# logging.error("This is a error log")
# logging.critical("This is a critical log")

class Raw_Reading_WithoutProcessBar:
    def __init__(self,file_psm1,file_psm2,Save_path,raw1,raw2):

        rawfile1 = MSFileReader(raw1,controllerType='MS',massAnalyzerType='FTMS',activationType='HCD',detectorType='HCD')
        rawfile2 = MSFileReader(raw2,controllerType='MS',massAnalyzerType='FTMS',activationType='HCD',detectorType='HCD')
        d1 = Path(file_psm1)
        d2 = Path(file_psm2)
        Save = Path(Save_path)
        fo = []
        d = {1: d1, 2: d2}
        Intensity_final_use = {}
        Peptide_final_use = {}
        AnalysisSignal.write(f"Start Processing Raw LC-MS Data\nRAW1: {raw1}\nRAW2: {raw2}")
        logging.info(f"Start Processing Raw LC-MS Data\nRAW1: {raw1}\nRAW2: {raw2}")
        for i, f in d.items():
            # print(f)
            file_path = f
            print(file_path)
            AnalysisSignal.write(f"xlsx/tsv File:{file_path}")
            logging.info(f"xlsx/tsv File:{file_path}")
            # Get the parent directories of each file
            parent_directories = (os.path.dirname(file_path))
            # Remove characters before and including the last part of each directory
            modified_directories = parent_directories.split("\\")[-1]
            modified_dir = modified_directories

            form = pd.read_csv(file_path, sep='\t', header=0)
                # watchpoints.watch(counter)
            folder = os.path.join(Save, f"{modified_dir}_RawDataPlot_{i}").replace("/", "\\")
            AnalysisSignal.write(folder)
            os.makedirs(folder)
            fo.append(folder)
            # 一次性进行数据筛选
            filtered_data = form[(form['Intensity'] != 0) & (form['MSFragger Localization'].notna())]
            Intensity = filtered_data['Intensity'].reset_index(drop=True)
            Intensity_final_use[i] = Intensity
            Peptide_final_use[i] = filtered_data['Peptide'].reset_index(drop=True)
        AnalysisSignal.write('\n')
        AnalysisSignal.write(Intensity_final_use)
        AnalysisSignal.write('\n')
        AnalysisSignal.write(Peptide_final_use)
        for p,Intensity in Intensity_final_use.items():
            # watchpoints.watch(p)
            if p == 1:
                rawfile = rawfile1
            else:
                rawfile = rawfile2
            # 获取第一张scan number
            scan_num_first = rawfile.GetFirstSpectrumNumber()
            # 获取最后一张scan number
            scan_num_last = rawfile.GetLastSpectrumNumber()
            dic = []
            Scan = []
            # 处理 scans
            for scan_num in range(scan_num_first, scan_num_last + 1):
                # AnalysisSignal.write(f'\nStarting Analysis Scan Number:{scan_num}')
                mass_list = rawfile.GetMassListFromScanNum(scan_num, scanFilter='')[0]
                mass = mass_list[0]
                intensity = mass_list[1]
                RT = rawfile.RTFromScanNum(scan_num)
                AnalysisSignal.write(f'Starting Analysis with RT:{RT} min')
                # logging.info(f'Starting Analysis with RT:{RT} min')
                for index1, a in enumerate(Intensity):
                    for index2, item in enumerate(intensity):
                        if item == a:
                            peptide = Peptide_final_use[p][index1]
                            Mass_F = mass[index2]
                            Scan.append(scan_num)
                            # watchpoints.watch(Scan)
                            dic.append((RT, Mass_F, peptide, item))
                            AnalysisSignal.write(f"RT:{RT}\nPeptide:{peptide}\nMass:{Mass_F}")
                            logging.info(f"RT:{RT}\nPeptide:{peptide}\nMass:{Mass_F}")
            AnalysisSignal.write("First Part Done! Start the second Part")
            logging.info("First Part Done! Start the second Part")
            # 处理结果
            # print(dic)
            # print(Scan)
            # 设置默认字体为 Times New Roman
            plt.rcParams['font.family'] = 'Times New Roman'
            # 处理 images
            for num in Scan:
                # watchpoints.watch(num)
                RTf = rawfile.RTFromScanNum(num)
                AnalysisSignal.write(f'\nStarting Ploting with RT:{RTf}')
                RTf_formatted = f"{RTf:.3f}"
                # 设置图形大小为 1800x1100 像素
                plt.figure(figsize=(18, 11))
                mass_list = rawfile.GetMassListFromScanNum(num, scanFilter='')[0]

                plt.plot(mass_list[0], mass_list[1], linestyle='-')  # 去掉 marker 参数
                # 记录每个 mass 的最高强度和对应的文本标签
                mass_labels = {}
                for m, mass_value in enumerate(mass_list[0]):
                    # watchpoints.watch(mass_value)
                    intensity_value = mass_list[1][m]
                    for item in dic:
                        if mass_value == item[1]:
                            peptide_label = f"Peptide: {item[2]}\nMass:{item[1]:.2f}\nIntensity: {intensity_value:.2f}"
                            if mass_value not in mass_labels or intensity_value > mass_labels[mass_value][1]:
                                mass_labels[mass_value] = (peptide_label, intensity_value)
                # 添加文本标签
                for mass_value, (peptide_label, intensity_value) in mass_labels.items():
                    plt.text(mass_value, intensity_value, peptide_label, fontsize=8, ha='center', va='bottom')
                plt.title(f"Retention Time: {RTf_formatted}")
                plt.xlabel("Mass (m/z)")
                plt.ylabel("Intensity")

                # 设置 Y 轴限制从 0 开始
                plt.ylim(0, )  # 根据需要调整 ymax
                plt.savefig(f"{fo[p-1]}/Scan_{num}.png", dpi=1000)
                plt.close()
        AnalysisSignal.write("The Plot Part of Analysis Thread finished!")
        logging.info("The Plot Part of Analysis Thread finished!")