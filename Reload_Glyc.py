import threading
import pandas as pd
import os
import re
from AnalysisSignalModule import AnalysisSignal
import logging
outputSignal = AnalysisSignal()
class GDASConnectionFasta_Reload(threading.Thread):
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

    def __init__(self,psm1_path,psm2_path,fasta_file):
        super().__init__()
        self.total_output = []
        self.d = [psm1_path, psm2_path]
        self.both_ID = []
        self.fasta_file = fasta_file

    def run(self):
        for p,psm_path in enumerate(self.d):
            if type(psm_path) == str:
                path = psm_path
                outputSignal.write(path)
                # self.update_text_browser1(f"{path}")
                parent_directories = (os.path.dirname(path))
                # Remove characters before and including the last part of each directory
                modified_directories = parent_directories.split("\\")[-1]
                modified_dir = modified_directories.split("/")[-1]
                filename = modified_dir

                Total1 = []

                i = os.path.join(path)
                print(i)
                outputSignal.write(i)
                form = pd.read_csv(i)
                # 筛选数据
                filtered_data = form['protein_name']
                # 检查筛选后的数据是否为空
                if not filtered_data.empty:
                    # 获取筛选后的数据中的"Protein ID"列的值
                    protein_ID = filtered_data.tolist()
                    print(protein_ID)
                    outputSignal.write(protein_ID)
                else:
                    print("No rows found where 'protein_name' is not empty.")
                    outputSignal.write("No rows found where 'protein_name' is not empty.")
                print(form)
                outputSignal.write(form)
                list1 = protein_ID
                list2 = []
                for element in list1:
                    match = re.findall(r'\|([A-Z0-9]+)\|', element)
                    list2.append(match[0])
                list2 = set(list2)
                for item in list2:
                    Total1.append(item)
                Total2 = set(Total1)
                Total2 = list(Total2)
                # print(Total2)
                # outputSignal.write(Total2)
                for i in Total2:
                    self.both_ID.append(i)
                # outputSignal.write(Total2)
            elif type(psm_path) == list:
                for path in psm_path:
                    outputSignal.write(path)
                    Total1 = []
                    i = os.path.join(path)
                    print(i)
                    outputSignal.write(i)
                    form = pd.read_csv(i)
                    # 筛选数据
                    filtered_data = form['protein_name']
                    # 检查筛选后的数据是否为空
                    if not filtered_data.empty:
                        # 获取筛选后的数据中的"Protein ID"列的值
                        protein_ID = filtered_data.tolist()
                        print(protein_ID)
                        outputSignal.write(protein_ID)
                    else:
                        print("No rows found where 'protein_name' is not empty.")
                        outputSignal.write("No rows found where 'protein_name' is not empty.")
                    print(form)
                    outputSignal.write(form)
                    list1 = protein_ID
                    list2 = []
                    for element in list1:
                        match = re.findall(r'\|([A-Z0-9]+)\|', element)
                        list2.append(match[0])
                    list2 = set(list2)
                    for item in list2:
                        Total1.append(item)
                    Total2 = set(Total1)
                    Total2 = list(Total2)
                    # print(Total2)
                    # outputSignal.write(Total2)
                    for i in Total2:
                        self.both_ID.append(i)
                        # outputSignal.write(Total2)
        del Total2,Total1,list1,list2
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
                            sy = GDASConnectionFasta_Reload.find_index_with_regex(fasta, item3)
                            output = ">" + fasta[sy]
                            print(output)
                            # r.write(output)
                            self.total_output.append(output)
                            outputSignal.write(output)

            # r.close()
        o.close()
        return self.total_output
        print("Done!")
        outputSignal.write('Done!')
        outputSignal.write('Fasta File Reload Done!')


    # 下面内容为修改已有的Fasta文件
    # for each in list_analysis:
    #     if each != 0:
    #         for i in range(len(list1)):
    #             list2.append(list1[i])
    #     else:
    #         print(each)
    # from requests.adapters import HTTPAdapter
    # from requests.packages.urllib3.util.ssl_ import create_urllib3_context
    # ORIGIN_CIPHERS = ('ECDH+AESGCM:DH+AESGCM:ECDH+AES256:DH+AES256:ECDH+AES128:DH+AES:ECDH+HIGH:'
    #                   'DH+HIGH:ECDH+3DES:DH+3DES:RSA+AESGCM:RSA+AES:RSA+HIGH:RSA+3DES')
    # class DESAdapter(HTTPAdapter):
    #     def __init__(self, *args, **kwargs):
    #         """
    #         A TransportAdapter that re-enables 3DES support in Requests.
    #         """
    #         import random
    #         CIPHERS = ORIGIN_CIPHERS.split(':')
    #         random.shuffle(CIPHERS)
    #         CIPHERS = ':'.join(CIPHERS)
    #         self.CIPHERS = CIPHERS + ':!aNULL:!eNULL:!MD5'
    #         super().__init__(*args, **kwargs)
    #
    #     def init_poolmanager(self, *args, **kwargs):
    #         context = create_urllib3_context(ciphers=self.CIPHERS)
    #         kwargs['ssl_context'] = context
    #         return super(DESAdapter, self).init_poolmanager(*args, **kwargs)
    #
    #     def proxy_manager_for(self, *args, **kwargs):
    #         context = create_urllib3_context(ciphers=self.CIPHERS)
    #         kwargs['ssl_context'] = context
    #         return super(DESAdapter, self).proxy_manager_for(*args, **kwargs)
    # import bs4
    # import requests.adapters
    # import logging
    # import requests
    # # Defination a Function: 通过uniprotID获取页面内容
    # def get_result(ID):
    #     headers = {
    #         'User-Agent': 'User-Agent:Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/56.0.2924.87 Safari/537.36'}
    #     # 解决报错requests.exceptions.SSLError: HTTPSConnectionPool
    #     logging.captureWarnings(True)
    #     s = requests.Session()
    #     host = "https://www.uniprot.org/uniprotkb/"
    #     url = host + ID + '.fasta'
    #     s.mount(host, adapter=DESAdapter())
    #     res = s.get(url, headers=headers)
    #
    #     # host = "https://www.uniprot.org/uniprot/"
    #     # url = host+ID
    #     # res = requests.get(url, headers=headers)
    #
    #     soup = bs4.BeautifulSoup(res.text, "html.parser")
    #     soup.get_text()
    #     result = str(soup)
    #     return result
    # with open("OPE_test_8ul_UniProt_Download.fasta",'w',encoding='utf-8') as f:
    #     # for item in Total2:
    #     modification = '>'
    #     a = get_result(item)
    #     t = a.replace('&gt;',modification)
    #     print(t)
    #     f.write(t)
    # f.close()