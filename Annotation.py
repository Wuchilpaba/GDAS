# import win32com.client
import pandas as pd
import re,os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from pyteomics import mzml
from pymsfilereader import MSFileReader
import pathlib
mpl.rcParams['font.family'] = 'Times New Roman'
mapping = {
    'NeuAc': 'S',
    'Hex': 'H',
    'HexNAc': 'N',
    'Fuc': 'F',
    'NeuGc': 'G'
}
def check_h_number(glycan_string):
    """Judge the High-mannose or not"""
    match = re.search(r'H(\d+)', glycan_string)
    if match:
        h_number = int(match.group(1))
        return h_number >= 5
    return False
def convert_glycan_to_abbreviation(single_glycan, mapping):
    """
    将单个聚糖结构转换为缩写形式

    Parameters:
    single_glycan (str): 单个聚糖结构字符串
    mapping (dict): 糖组分映射字典

    Returns:
    str: 缩写后的聚糖结构
    """
    # 如果single_glycan为空，返回空字符串
    if pd.isna(single_glycan) or single_glycan == '':
        return ''

    # 提取各组分数量
    components = {}
    for glycan_type, abbrev in mapping.items():
        pattern = f'{glycan_type}\\((\\d+)\\)'
        matches = re.findall(pattern, single_glycan)
        if matches:
            components[abbrev] = int(matches[0])

    # 按照S-H-N-F的顺序构建缩写字符串
    ordered_components = ['S', 'G', 'H', 'N', 'F']
    result_parts = []

    for comp in ordered_components:
        if comp in components and components[comp] > 0:
            result_parts.append(f'{comp}{components[comp]}')

    return ''.join(result_parts)

def GDAS_result(path):
    protein = pd.read_excel(path,sheet_name='Sheet1')['Protein'].tolist()
    return protein

def clean_column_names(df):
    """清理列名中的换行符和多余空格"""
    df.columns = [re.sub(r'\s+', ' ', col.replace('\n', ' ').strip()) for col in df.columns]
    return df

def Byonic_result(path_list,protein_list):
    total = {}
    Glycan = {}
    Peptide = {}
    for index,path in enumerate(path_list):
        df = pd.read_excel(path, sheet_name='Spectra', engine="openpyxl")
        # 清理列名中的换行符和多余空格
        df = clean_column_names(df)
        try:
            mask = df['Protein Name'].apply(
                lambda x: any(protein in str(x) for protein in protein_list)
            )
            df_clean = df[mask].rename(columns={
                'Peptide < ProteinMetrics Confidential >': 'Peptide',
                'Glycans NHFAGNa': 'Glycans',
                'Scan #': 'Scan'
            })
            protein_data = df_clean
        except:
            protein_data = pd.DataFrame
        for protein in protein_list:
            total[protein] = {}
            Glycan[protein] = {}
            Peptide[protein] = {}
            total[protein][index] = []
            Glycan[protein][index] = []
            Peptide[protein][index] = []
            filtered_data0 = protein_data[protein_data['Protein Name'].str.contains(protein, na=False)]
            filtered_data1 = filtered_data0.sort_values(by=['Score'], ascending=False).reset_index(drop=True)
            if not filtered_data1.empty:
                scan_list = filtered_data1['Scan'].tolist()
                peptide_list = filtered_data1['Peptide'].tolist()
                glycan_list = filtered_data1['Glycans'].tolist()
            else:
                raise ValueError('Fatal Error: During Ploting Annotation, result in Byonic data pre-process')
            n_glycan_list = []
            for glycan in glycan_list:
                n_glycan_list.append(convert_glycan_to_abbreviation(glycan,mapping))
            glycan_list = n_glycan_list
            count = 0
            i = 0
            while count < 10 and i < len(scan_list):
                scan = scan_list[i]
                p = peptide_list[i]
                g = glycan_list[i]

                if g != '':
                    secondary = scan.split(' ')
                    for se in secondary:
                        if 'scan' in se:
                            s = se.split('=')
                            total[protein][index].append(s[1])
                    Glycan[protein][index].append(g)
                    Peptide[protein][index].append(p)
                    count += 1
                i += 1
    return total,Peptide,Glycan

def plot(mz,intensity,scan,rt,peptide,glycan,protein,path,filename):

    a4_width, a4_height = 8.27, 11.69
    fig = plt.figure(figsize=(a4_width, a4_height))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.5])
    # 文字区
    ax_text = fig.add_subplot(gs[0])
    ax_text.axis("off")

    text_info = (
        "MS Annotation\n"
        f"File: {filename}\n"
        f"Scan: {scan}\n"
        f"Scan Time: {rt}min\n"
        f"Number of peaks: {len(mz)}\n"
        f"Protein: {protein}\n"
        f"Peptide: {peptide}\n"
        f"Glycan: {glycan}"
    )
    ax_text.text(0, 1, text_info, va="top", ha="left", fontsize=12,transform=ax_text.transAxes)

    # 质谱图
    ax_plot = fig.add_subplot(gs[1])
    ax_plot.vlines(mz, [0], intensity, color="black", linewidth=0.8)
    ax_plot.set_xlabel("m/z")
    ax_plot.set_ylabel("Intensity")
    ax_plot.set_title(f"MS Spectrum (Scan {scan}, RT={rt:.2f} min)")
    ax_plot.spines['top'].set_visible(False)
    ax_plot.spines['right'].set_visible(False)
    ax_plot.set_xlim(left=0)
    ax_plot.set_ylim(bottom=0)
    ax_plot.ticklabel_format(style="sci", axis="y", scilimits=(0, 0),useMathText=True)
    plt.tight_layout()
    originalPath = pathlib.Path(path)
    plt.savefig(originalPath / f"{filename}_{scan}.png", dpi=600)
    # plt.show()

def main(impact_path,path_list,path_save_plot,MSFileList):
    '''impact_path: GDAS result excel file path, path_list: Byonic file path list, path_save_plot: output file path'''
    predict = GDAS_result(impact_path)
    Total,Peptide,Glycan = Byonic_result(path_list,predict)
    for MSFile in MSFileList:
        for protein in Total:
            ext = os.path.splitext(MSFile)[1].lower()
            filename = os.path.splitext(os.path.basename(MSFile))[0]
            if ext == ".raw":
                ms_file = MSFileReader(filename=MSFile)
                # ms_file = win32com.client.Dispatch("MSFileReader.XRawfile.1")
                # ms_file.Open(MSFile)
                for i,scan in Total[protein].items(): # 第i个文件
                    # for s in scan:
                    gly = Glycan[protein][i]
                    for n,can in enumerate(gly):
                        s = int(scan[n])
                        if 'F' in can or 'S' in can or check_h_number(can):
                            rt = ms_file.RTFromScanNum(s)
                            mz,intensity = ms_file.GetMassListFromScanNum(s)[0]
                            # mz = tu[0]
                            # intensity = tu[1]
                            plot(mz, intensity, s, rt, Peptide[protein][i][n], can, protein, path_save_plot, filename)
                        else:
                            pass
            elif ext == ".mzml":
                reader = mzml.MzML(MSFile)
                spec = None
                for i,scan in Total[protein].items():
                    gly = Glycan[protein][i]
                    for n,can in enumerate(gly):
                        s = int(scan[n])
                        if 'F' in can or 'S' in can or check_h_number(can):
                            for spectrum in reader:
                                sid = spectrum.get("id", "")
                                if f"scan={s}" in sid:
                                    spec = spectrum
                                    break
                            if spec is None:
                                raise ValueError(f"In {MSFile}, not found ScanNum={s}")
                            mz = spec["m/z array"]
                            intensity = spec["intensity array"]
                            rt = spec["scanList"]["scan"][0]["scan start time"]
                            plot(mz, intensity, s, rt, Peptide[protein][i][n], can, protein, path_save_plot,filename)
                        else:
                            pass