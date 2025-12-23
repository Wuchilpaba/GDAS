import modin.pandas as pd
import pandas
import ray
import logging
import traceback
import pickle
import os
import re


def parse_index_range(index_range_str):
    indices = re.findall(r'\d+', index_range_str)
    return int(indices[0]), int(indices[1])

def find_o_glycosylation_sites(seq, index_range_strings):
    if seq and index_range_strings:
        start_index, end_index = parse_index_range(index_range_strings)
        o_glyco_sites = [
            (m.group()[0], m.start() - sum(len(ann) for ann in re.findall(r'\[.*?\]', seq[:m.start()])) + start_index)
            for m in re.finditer(r'[ST]\[O-Glycosylation.*?\]', seq)
        ]
        cleaned_seq = re.sub(r'\[.*?\]', '', seq)
        return o_glyco_sites, cleaned_seq
    return None
def is_close_float(a, b, tolerance=1e-6):
    """
    Compare two floating-point numbers to see if they are equal within a given margin of error.
    :param a: First Float Number
    :param b: Second Float Number
    :param tolerance: tolerance (allowed error)
    :return: Boolean value indicating whether two numbers are close
    """
    return abs(a - b) <= tolerance


@ray.remote
def find_most_similar_sequence(PSM_DataFrame, Peak_DataFrame):

    try:
        # 提取所需的列
        possiblechange_columns = ['Scan Retention Time', 'Retention Time']
        valid_column = None
        for col in possiblechange_columns:
            if col in PSM_DataFrame.columns:
                valid_column = col
                break
        peak_data = Peak_DataFrame[['Full Sequence', 'MS2 Retention Time',
                                    'Peptide Monoisotopic Mass', 'Protein Group', 'Peak intensity']]
        psm_data = PSM_DataFrame[['Full Sequence', valid_column,
                                  'Peptide Monoisotopic Mass', 'Protein Accession',
                                  'GlycanMass', 'Plausible GlycanComposition', 'Start and End Residues In Protein']]

        # 创建匹配字典
        psm_dict = {
            (row['Full Sequence'], row['Protein Accession']): row
            for _, row in psm_data.iterrows()
        }

        final_results = []
        for _, peak_row in peak_data.iterrows():
            key = (peak_row['Full Sequence'], peak_row['Protein Group'])

            if key in psm_dict:
                psm_row = psm_dict[key]
                # 检查浮点数匹配
                if is_close_float(peak_row['MS2 Retention Time'], psm_row[valid_column]) and \
                    is_close_float(peak_row['Peptide Monoisotopic Mass'], psm_row['Peptide Monoisotopic Mass']):
                    site, cleaned_seq = find_o_glycosylation_sites(peak_row['Full Sequence'],psm_row['Start and End Residues In Protein'])
                    final_results.append({
                        'ID': peak_row['Protein Group'],
                        'Peptide': cleaned_seq,
                        'Peptide with Site': peak_row['Full Sequence'],
                        'Glycan': psm_row['Plausible GlycanComposition'],
                        'GlycanMass': float(psm_row['GlycanMass']),
                        'Intensity': float(peak_row['Peak intensity']),
                        'Site': site
                    })

        return final_results

    except Exception as e:
        logging.error(f"Error in remote task: {e}")
        traceback.print_exc()
        return []


if __name__ == "__main__":
    try:
        current_directory = os.getcwd()
        worker_ray_path = os.path.join(current_directory)
        ray.init(ignore_reinit_error=True)

        with open(rf'{worker_ray_path}/data_cache_Meta_Intensity.pkl', 'rb') as f:
            data = pickle.load(f)
        PSM_dataframe, Peak_dataframe = data
        PSM_dataframe = PSM_dataframe.applymap(lambda x: x.strip() if isinstance(x, str) else x)
        Peak_dataframe = Peak_dataframe.applymap(lambda x: x.strip() if isinstance(x, str) else x)
        # Version Possible change in result files
        possiblechange_columns = ['Scan Retention Time', 'Retention Time']
        valid_column = None
        for col in possiblechange_columns:
            if col in PSM_dataframe.columns:
                valid_column = col
                break

        # 转换数据类型，避免字符串类型参与浮点数运算
        PSM_dataframe[valid_column] = pd.to_numeric(PSM_dataframe[valid_column], errors='coerce')
        PSM_dataframe['Peptide Monoisotopic Mass'] = pd.to_numeric(PSM_dataframe['Peptide Monoisotopic Mass'],
                                                                   errors='coerce')

        Peak_dataframe['MS2 Retention Time'] = pd.to_numeric(Peak_dataframe['MS2 Retention Time'], errors='coerce')
        Peak_dataframe['Peptide Monoisotopic Mass'] = pd.to_numeric(Peak_dataframe['Peptide Monoisotopic Mass'],
                                                                    errors='coerce')
        # 并行任务
        chunk_size = len(Peak_dataframe) // 4 + 1
        tasks = [
            find_most_similar_sequence.remote(PSM_dataframe, Peak_dataframe.iloc[i:i + chunk_size])
            for i in range(0, len(Peak_dataframe), chunk_size)
        ]
        results = ray.get(tasks)
        final_results = [item for sublist in results for item in sublist]
        print('Final Results:', final_results)
        final_df = pandas.DataFrame(final_results)
        current_directory = os.getcwd()
        cache_file_path = os.path.join(current_directory, 'cache_file_task_Meta.pkl')
        final_df.to_pickle(cache_file_path)
        print(f"Final results have been saved to {cache_file_path}")
    except Exception as e:
        logging.error(f"Error during multiprocessing: {e}")
        traceback.print_exc()
    finally:
        print('Finished MetaMorpheus Intensity pairing processing with Ray')
