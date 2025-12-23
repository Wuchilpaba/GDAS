import ctypes
import json
import os, traceback, subprocess, logging
import pathlib
import re
import time,pickle,sys, shutil
base_path = os.path.abspath(os.path.dirname(__file__))
internal_path = os.path.join(base_path, '_internal', 'psutil')
if internal_path not in sys.path:
    sys.path.insert(0, internal_path)
try:
    import psutil
    print(f"psutil version: {psutil.__version__}")
    print(f"Loaded from: {psutil.__file__}")
except ImportError as e:
    print(f"Failed to load psutil: {e}")
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QTextCursor
from PyQt5.QtCore import QObject, pyqtSignal, QEventLoop, QTimer, QThread, pyqtSlot, QRunnable, QThreadPool
from PyQt5.QtWidgets import QInputDialog, QMessageBox, QGraphicsScene, QGraphicsPixmapItem
from math import log10,log2
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QFileDialog, QVBoxLayout, QWidget, QToolTip, QMenu
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QCursor

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from MSFragger_Fasta_Reload_Script import GDASConnectionFasta_Reload as Reload_MSFragger
from MSFragger_Fasta_Reload_Script import outputSignal as MSFragger_reload_signal
from MSFragger_Analysis import Analysis as MSFragger_Process
from MSFragger_Analysis import analysis_signal as MSFragger_raw_reading_signal

from GDAS_GUI import Ui_MainWindow

from Analysis_GlycReSoft import Analysis_GlycReSoft as glyc_analysis_process
from Analysis_GlycReSoft import analysis_signal as GlyAnalysisoutput

from MataMorpheus_Reload import MetaMorpheus_Fasta_Reload
from MetaMorpheus_Analysis import Analysis_MetaMorpheus as MetaMorpheusAnalysisProcess
from MetaMorpheus_Analysis import MetaMorpheus_output as MetaMorpheusAnalysisoutput
from MataMorpheus_Reload import outputSignal as MetaReloadOutput

from Reload_Glyc import GDASConnectionFasta_Reload as GlyReSoft
from Reload_Glyc import outputSignal as GlyreloadOutput

from data_module import DataProvider

import docker

from PCA_Glycosylation_Analysis import final_signal
import PCA_Glycosylation_Analysis
from resource_path import resource_path

# Preparing the Program
myappid = "myApplication"
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
with open("main", mode='w') as f:
    pass
f.close()
# logger = Logger.Logger(level='debug').logger
# logger.info("Logger initialized successfully.")
CONFIG_FILE = resource_path('resources/config.json')
script_root_path = pathlib.Path(__file__).resolve().parent
All_Ca = script_root_path / "All_Cache"
if not All_Ca.exists():
    All_Ca.mkdir(parents=True, exist_ok=True)
else:
    for item in All_Ca.iterdir():
        if item.is_file() or item.is_symlink():
            item.unlink()  # Delete files or symbolic links
            logging.info(f"Deleted: {item}")
        elif item.is_dir():
            shutil.rmtree(item)  # Delete subfolders recursively
            logging.info(f"Deleted: {item}")
def load_config(config_file):
    if os.path.exists(config_file):
        with open(config_file, 'r') as file:
            return json.load(file)
    return {}
def save_config(config_data):
    with open(CONFIG_FILE, 'w') as file:
        json.dump(config_data, file, indent=4)

class Signal(QObject):
    text_update = pyqtSignal(str)
    def write(self, text):
        self.text_update.emit(str(text))
        loop = QEventLoop()
        QTimer.singleShot(100, loop.quit)
        loop.exec_()
        QApplication.processEvents()
class MatplotlibCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig, self.ax = plt.subplots()
        super(MatplotlibCanvas, self).__init__(fig)
        self.setParent(parent)
        self.data_points = []

    def plot_volcano(self, log2fc, p_values, labels, fcvalue, pvalue):
        self.data_points = list(zip(log2fc, p_values, labels))

        self.ax.clear()
        colors = []
        for log2fc_val, p_value in zip(log2fc, p_values):
            if log2fc_val > log2(float(fcvalue)) and p_value > -log10(float(pvalue)):
                colors.append('#FF0000')  # Red
            elif log2fc_val < -log2(float(fcvalue)) and p_value > -log10(float(pvalue)):
                colors.append('#0000FF')  # Blue
            else:
                colors.append('#A9A9A9') # Grey
        self.ax.scatter(log2fc, p_values, c=colors)

        # 标记阈值线，例如 log2FC > 1 或 log2FC < -1，p-value < 0.05
        self.ax.axhline(y=-log10(float(pvalue)), color='red', linestyle='--')
        self.ax.axvline(x=log2(float(fcvalue)), color='red', linestyle='--')
        self.ax.axvline(x=-log2(float(fcvalue)), color='red', linestyle='--')

        self.ax.set_xlabel('log2(Fold Change)')
        self.ax.set_ylabel('-log10(P-value)')

        max_abs_log2fc = max(abs(fc) for fc in log2fc)
        self.ax.set_xlim(-max_abs_log2fc - 1, max_abs_log2fc + 1)

        # 鼠标悬停事件
        self.mpl_connect('motion_notify_event', self.on_hover)
        self.draw()

    def on_hover(self, event):
        if event.inaxes == self.ax:
            for point in self.data_points:
                log2fc, p_value, label = point
                x = log2fc
                y = p_value
                if abs(event.xdata - x) < 0.1 and abs(event.ydata - y) < 0.1:
                    QToolTip.showText(QCursor.pos(), f"ID: {label}\nLog2(FC): {log2fc:.2f}\n-Log10(P-value): {p_value:.2e}")
                    break
            else:
                QToolTip.hideText()

    def save_figure(self, filename):
        self.figure.savefig(filename)
class ImageWindow(QMainWindow):
    def __init__(self, log2fc, p_values, labels, fcvalue, pvalue):
        super(ImageWindow, self).__init__()

        self.canvas = MatplotlibCanvas()
        self.init_ui()
        self.canvas.plot_volcano(log2fc, p_values, labels, fcvalue, pvalue)

    def init_ui(self):
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(resource_path("resources/CCMS.tif")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        layout.addWidget(self.canvas)
        self.setCentralWidget(central_widget)

        self.setWindowTitle('Volcano Plot Viewer')
        self.resize(800, 600)

        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)

    def show_context_menu(self, pos):
        context_menu = QMenu(self)
        save_as_action = QAction('Save As', self)
        save_as_action.triggered.connect(self.save_as)
        context_menu.addAction(save_as_action)
        context_menu.exec_(self.mapToGlobal(pos))
    def save_as(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getSaveFileName(self, 'Save Image As', '',
                                                   'PNG Files (*.png);;JPEG Files (*.jpg);;All Files (*)')
        if file_path:
            self.canvas.save_figure(file_path)
class MetaProcessingTask(QRunnable):
    def __init__(self, MS_File, result, binds, Setting, signal_emitter):
        super().__init__()
        self.MS_File = MS_File
        self.result = result
        self.binds = binds
        self.Setting = Setting
        self.emitter = signal_emitter
    def run(self):
        try:
            # Start_Glyc_container function
            infor = self.start_opair_container(self.MS_File,self.result,self.binds, self.Setting)
            # wait for container
            co = infor.get('exit_code')
            logs = infor.get('logs')
            # AnalysisThread.update_signal[str].emit(): argument 1 has unexpected type 'dict'
            d = [str(co),str(logs)]
            update = ",".join(d)
            self.emitter.emit(update)
        except Exception as e:
            self.emitter.emit(str(e))
    def start_opair_container(self, MS_File, result, binds, set):
        '''
        :param MS_File: Design to fill with MS files, format: mzml or raw, just for name of files
        :param result: Result folder path
        :param binds: file paths which needs to be loaded in docker container
        :param set: MataMorpheus Setting Folder
        :return: None, this function is designed to work and that means it would not return
        '''

        print('Starting O-Pair processing in Docker')
        MS_File = pathlib.Path(MS_File)
        client = docker.from_env()

        all_volumes = {
            binds: {  # fasta+Rules+MSFiles
                "bind": "/All_Cache",
                "mode": "rw"
            },
            set: {
                "bind": "/MetaMorpheusAdditionalSettings",
                "mode": "rw"
            },
            result: {  # workdir
                "bind": "/mnt",
                "mode": "rw"
            }
        }

        name = MS_File.name
        container = client.containers.run(
            image="smithchemwisc/metamorpheus:latest",
            detach=True,
            tty=True,
            auto_remove=True,
            volumes=all_volumes,
            command=[
                '-t', '/All_Cache/Cache.toml',
                '-s', f'/All_Cache/{name}',
                '-d', '/All_Cache/database_O-Pair_Tem.fasta',
                '-o', '/mnt',
                '--mmsettings', '/MetaMorpheusAdditionalSettings'
            ],
            name=f"{name}_O-Pair_Docker"
        )
        try:
            exit_code = container.wait()
            logs = container.logs().decode('utf-8')
            return {
                'exit_code': exit_code,
                'logs':logs,
                'container_id':container.id
            }
        finally:
            try:
                container.remove()
            except:
                pass
class GlycProcessingTask(QRunnable):
    def __init__(self,MS_File, result, binds, signal_emitter):
        super().__init__()
        self.MS_File = MS_File
        self.result = result
        self.binds = binds
        self.emitter = signal_emitter
    def run(self):
        try:
            # Start_Glyc_container function
            infor = self.start_glyc_container(self.MS_File,self.result,self.binds)
            # wait for container
            exit_status = infor.get('exit_code')
            logs = infor.get('logs')
            d = [str(exit_status), str(logs)]
            update = ",".join(d)
            self.emitter.emit(update)
        except Exception as e:
            self.emitter.emit(str(e))
    def start_glyc_container(self, MS_File, result, binds):
        '''
        :param MS_File: Design to fill with MS files, format: mzml. GlycReSoft requires the format must be mzml
        :param result: Result folder path
        :param binds: Other file which needs to be loaded in docker container
        :return: None, this function is designed to work and that means it would not return
        '''
        print('Starting GlycReSoft processing in Docker')

        MS_File = pathlib.Path(MS_File)
        client = docker.from_env()

        all_volumes = {
            binds: {  # fasta+Rules+MSFiles
                "bind": "/All_Cache",
                "mode": "rw"
            },
            result: {  # workdir
                "bind": "/mnt",
                "mode": "rw"
            }
        }

        name = MS_File.name
        container = client.containers.run(
            image="mobiusklein/glycresoft:v0.4.24",
            detach=True,
            tty=True,
            auto_remove=True,
            volumes=all_volumes,
            command=[
                '/bin/bash', '-c',
                f'glycresoft build-hypothesis glycan-combinatorial /All_Cache/rules_file.txt combinatorial-database -n "Combinatorial Human N-Glycans" && '
                f'glycresoft build-hypothesis glycopeptide-fa -g /All_Cache/rules_file.txt -s combinatorial -e trypsin -c "Carbamidomethyl (C)" -v "Oxidation (M)" -v "Pyro-glu from Q (Q@N-term)" -v "Pyro-glu from E (E@N-term)" -n "Cache" "/All_Cache/database_Glyc_Tem.fasta" "fasta-glycopeptides.db" && '
                f'glycresoft mzml preprocess -an glycopeptide -tn 10 -mn 1 -c 12 -n {name} "/All_Cache/{name}" "/mnt/{name}" && '
                f'glycresoft analyze search-glycopeptide -n {name} --export csv fasta-glycopeptides.db /mnt/{name} 1 -o "/mnt/Cache-glycopepitides-in-samples_{name}.db" && '
                'echo \'Commands executed successfully\''
            ],
            name=f"{name}_Glyc_Docker"
        )
        try:
            exit_code = container.wait()
            logs = container.logs().decode('utf-8')
            return {
                'exit_code':exit_code,
                'logs': logs,
                'container_id':container.id
            }
        finally:
            try:
                container.remove()
            except:
                pass
class AnalysisThread(QThread):
    update_signal = pyqtSignal(str)
    volcano_data_ready = pyqtSignal(list, list, list, float, float)
    debugging_signal = pyqtSignal(str)
    def __init__(self, file_raw1, file_raw2, file_psm1, file_psm2, save_path, output_signal, analysis_signal, reload_signal, fasta1_filepath, choice,value_2, FcValue, PValue, ContainerThread,Annotation):
        super().__init__()
        self.stopped = False
        logging.info("Analysis Thread Started")
        self.file_raw1 = file_raw1
        self.file_raw2 = file_raw2
        self.file_psm1 = file_psm1
        self.file_psm2 = file_psm2
        self.save_path = save_path
        # self.raw1 = file_raw1
        # self.raw2 = file_raw2
        self.output_signal = output_signal
        self.analysis_signal = analysis_signal
        # self.AnalysisSignal = raw_reading_signal
        self.outputSignal = reload_signal
        self.fasta1_filepath = fasta1_filepath
        self.choice = choice

        self.value2 = value_2 # Mode in MSF
        self.FcValue = FcValue # FC
        self.PValue = PValue
        self.containerThread = ContainerThread
        self.Annotation = Annotation
    # def start_MSF_container(self ,work_path ,result, manifest,fasta):
    #     import docker
    #     client = docker.from_env()
    #     print("Starting Ray container...")
    #     # '_internal'
    #     # Origin root path
    #     o = pathlib.Path.cwd()
    #     resources = pathlib.Path(o, 'resources')
    #     workflow = resources / 'cache.workflow'
    #     tools_path = resources / 'tools'
    #     # MSFile
    #     dynamic_binds = {}
    #     for path in work_path:
    #         path_obj = pathlib.Path(path) if not isinstance(path, pathlib.Path) else path
    #         if path_obj.exists() and path_obj.is_file():
    #             # 将文件挂载到容器目标路径
    #             dynamic_binds[str(path_obj.resolve())] = {
    #                 "bind": f"/home/{path_obj.name}",
    #                 "mode": "rw",
    #             }
    #     fixed_volumes = {
    #         manifest: {  # manifest
    #             "bind": "/home/cache.fp-manifest",
    #             "mode": "rw"
    #         },
    #         fasta: {    # fasta
    #             "bind": "/home/database_with_decoy.fasta",
    #             "mode": "rw"
    #         },
    #         result: {  # workdir
    #             "bind": "/mnt",
    #             "mode": "rw"
    #         },
    #         workflow: {  # workflow
    #             "bind": "/opt/cache.workflow",
    #             "mode": "rw"
    #         },
    #         tools_path: {
    #             "bind": "/tools",
    #             "mode": "ro"
    #         }
    #     }
    #     all_volumes = {**fixed_volumes, **dynamic_binds}
    #     # work_flow_path = pathlib.Path(inpath, '241211.workflow')
    #     container = client.containers.run(
    #         image="fcyucn/fragpipe:22.0",
    #         detach=True,
    #         tty=True,
    #         auto_remove=True,
    #         volumes=all_volumes,
    #         command=[
    #             "/bin/bash", "-c",
    #             "cd /fragpipe_bin/fragPipe-22.0/fragpipe/bin && "
    #             "./fragpipe --headless "
    #             "--workflow /opt/cache.workflow "
    #             "--manifest /home/cache.fp-manifest "
    #             "--workdir /mnt "
    #             "--config-tools-folder /tools && "
    #             "echo 'Commands executed successfully'"
    #         ],
    #         name="MSF_Docker",
    #
    #     )
    #     print(f"Started container:{container.id}")
    #     print('Waiting for the Docker return')
    #     # reader = client.containers.get(container.id).logs(stream=True, follow=True)
    #     # for line in reader:
    #     #     print(line.decode('utf-8'),end='',flush=True)
    #     return container





    def stop(self):
        self.stopped = False
        sys.exit()

    def run(self):
        while not self.isInterruptionRequested():
            sys.stdout = self.output_signal
            QApplication.processEvents()
            self.process_step()
            break

    def process_step(self):
        self.perform_task_based_on_choice()

    def perform_task_based_on_choice(self):
        # ["N-Glycopeptide", "O-Glycopeptide"]
        if self.choice == 1 and self.value2 == "N-Glycopeptide":
            self.handle_msr_fragger_task()
            self.handle_glycresoft_task()
        elif self.choice == 2 and self.value2 == "O-Glycopeptide":
            self.handle_msr_fragger_task()
            self.handle_metamorpheus_task()
        else:
            self.handle_final_analysis()

    def handle_msr_fragger_task(self):
        # self.debugging_signal.emit("Start dealing with results from MSF")
        # data = []
        # o = pathlib.Path.cwd()
        # inpath = pathlib.Path(o, 'resources')
        # with open(pathlib.Path(inpath,'cache.workflow'), 'r') as file:
        #     for line in file:
        #         line = line.strip()
        #         if line == '':
        #             # 保留空行
        #             data.append([line, None])  # 使用空字符串作为占位符
        #         elif line.startswith('#'):
        #             # 保留注释行，包括'#'符号
        #             data.append([line, None])  # 使用空字符串作为占位符
        #         else:
        #             parts = line.split('=', 1)
        #             if len(parts) == 2:
        #                 key, value = parts
        #                 # 将键和值添加到列表中的元组里
        #                 data.append([key.strip(), value.strip()])
        # df = pd.DataFrame(data, columns=['Key', 'Value'])
        # logging.info(f'Start handle Data with MSFragger in docker of fragpipe')
        # # 目标文件路径（原始Target序列）
        # target_fasta_path = pathlib.Path(inpath, 'database.fasta')
        # # 输出含Decoy的文件路径
        # decoy_fasta_path = pathlib.Path(inpath, 'database_with_decoy.fasta')
        # with open(self.fasta1_filepath, 'rb') as src_file:
        #     with open(pathlib.Path(inpath, 'database.fasta'), 'wb') as dst_file:
        #         dst_file.write(src_file.read())
        #     dst_file.close()
        # src_file.close()
        # with open(target_fasta_path, 'r') as source_file, open(decoy_fasta_path, 'w') as output_file:
        #     Fasta.write_decoy_db(source=source_file, output=output_file, file_mode='fasta')
        # output_file.close()
        # temp_output_path = decoy_fasta_path.with_suffix(".temp.fasta")
        # with open(decoy_fasta_path, 'r') as infile, open(temp_output_path, 'w') as outfile:
        #     for line in infile:
        #         if line.startswith('>'):
        #             line = line.replace("DECOY", "rev")
        #         outfile.write(line)
        # if decoy_fasta_path.exists():
        #     os.remove(decoy_fasta_path)
        # temp_output_path.rename(decoy_fasta_path)
        # # os.remove(temp_output_path)
        # print(f"New database with Decoy has been generated: {decoy_fasta_path}")
        # df.loc[df['Key'] == 'database.db-path', 'Value'] = '/home/database_with_decoy.fasta'
        # df = df.fillna(value='')
        # with open(pathlib.Path(inpath, 'cache.workflow'), 'w') as file:
        #     for index, row in df.iterrows():
        #         if pd.isnull(row['Key']) or row['Key'] == '':
        #             # 写入空行，写入了=None
        #             file.write("\n")
        #         elif row['Key'].startswith('#'):
        #             # 写入注释行
        #             file.write(f"{row['Key']}\n")
        #         else:
        #             # 写入键值对
        #             file.write(f"{row['Key']}={row['Value']}\n")
        # save_path = pathlib.Path(self.save_path) / 'MSFragger_Result'
        # save_path.mkdir(parents=True, exist_ok=True)
        # save = save_path
        # del save_path
        # # manifest file
        # manifestfile = inpath / 'cache.fp-manifest'
        # # fasta
        # fasta = inpath / 'database_with_decoy.fasta'
        # shutil.copy(manifestfile, save.parent)
        # shutil.copy(fasta, save.parent)

        # for index, each in enumerate(FileGroupCombine):
        #     new_records = []
        #     for indexi,file in enumerate(each):
        #         records = []
        #         file = pathlib.Path(file)
        #         name = file.name
        #         records.append(f"/home/{name}")
        #         if index == 0:
        #             records.append('Experiement')
        #         else:
        #             records.append('Control')
        #         records.append(str(indexi+1))
        #         records.append('DDA')
        #         new_records.append(records)
        #     try:
        #         with open(manifestfile, "w") as file:
        #             for record in new_records:
        #                 file.write("\t".join(record) + "\n")
        #         print(f"'{manifestfile}' updated successfully")
        #     except Exception as e:
        #         print(f"Error: {e}")
        #     # Key Process of MSF
        #     container = self.start_MSF_container(each,ResultCombine[index],manifestfile,fasta)
        #     exit_code = container.wait()
        #     logging.info(f"Container {container.id} finished with exit code {exit_code}")
        #     if exit_code == 0:
        #         logging.info("Command executed successfully, proceeding with next steps or files...")
        #         container.stop()
        #     else:
        #         logging.error("Command execution failed, check the container logs for details.")
        # After MSF
        # FileGroupCombine = [self.file_raw1, self.file_raw2]
        save1 = self.file_psm1
        save2 = self.file_psm2
        ResultCombine = [save1, save2]
        # pattern = re.compile(r'^psm\.tsv$')
        file_psm1_for_R_O = []
        file_psm2_for_R_O = []
        for i,Result in enumerate(ResultCombine):
            if isinstance(Result, list) and i == 0:
                for f in Result:
                    file_psm1_for_R_O.append(f)
            elif isinstance(Result, list) and i == 1:
                for f in Result:
                    file_psm2_for_R_O.append(f)
            else:
                raise TypeError('Input File Error from MSFragger Files')
        if file_psm1_for_R_O is [] or file_psm2_for_R_O is []:
            raise FileNotFoundError('After processing, result files return None or files could not be found')
        # MSFA = save / 'MSFragger_Analysis_Result'
        # Connect to the last analysis module of MSF
        save = pathlib.Path(self.save_path) / "MSFragger_Results_Analysis"
        self.MS_R_O = Reload_MSFragger.ret(file_psm1_for_R_O, file_psm2_for_R_O, self.fasta1_filepath, save, self.value2)
        MSProcess = MSFragger_Process(file_psm1_for_R_O, file_psm2_for_R_O, save, self.MS_R_O, self.FcValue, self.value2)
        Ratio3, self.MSF_fasta_output = MSProcess.run()
        logging.info(f'Result:{Ratio3}')
        self.update_signal.emit(f"{save}\\Protein_Ratio_Difference.png")
        i = 0
        j = 0
        for key in Ratio3.keys():
            if key[0] == 1:
                i += 1
            elif key[0] == 2:
                j += 1
        if i >= 3 and j >= 3:
            self.start_Statistic(Ratio3, self.save_path)
        else:
            logging.warning('Lack of sufficient data for model fitting')
        print('Analysis of MSF finished')

    def copy_file_to_folder(self, src_file: pathlib.Path, dst_folder: pathlib.Path):
        """Copy file from other path to Cache folder"""
        try:
            if src_file.exists():  # Check the resource
                if not dst_folder.exists():
                    dst_folder.mkdir(parents=True, exist_ok=True)
                dst_file = dst_folder / src_file.name  # Create the path of target file
                shutil.copy(src_file, dst_file)  # Copy
                print(f"{src_file} has copied to '{dst_file}'！")
            else:
                print(f"The file '{src_file}' does not exist！")
        except Exception as e:
            logging.error(f"failed to copy file {src_file}: {e}")
    def handle_glycresoft_task(self):
        print('Start handle data with GlycReSoft')
        logging.info(f'Start handle Data with GlycReSoft')
        script_root_path = pathlib.Path(__file__).resolve().parent
        All_Cache = script_root_path / "All_Cache"
        # self.debugging_signal.emit("Start dealing with results from Glyc")
        # Start Processing data in Docker
        P = [self.file_raw1, self.file_raw2]
        save = pathlib.Path(self.save_path) / 'GlycReSoft_Result'
        fasta = pathlib.Path(resource_path('resources')) / 'database_Glyc_Tem.fasta'
        rules = pathlib.Path(fasta).parent / "rules_file.txt"
        while not self.isInterruptionRequested():
            if len(self.MSF_fasta_output) == 0:
                self.debugging_signal.emit('Something unexpected has happened, please check the MS Files input. Task interrupted.')
                self.requestInterruption()
                break
            else:
                with open(fasta, 'w') as dst_file:
                    for sequence in self.MSF_fasta_output:
                        dst_file.write(sequence + "\n")
                break
        # dst_file.close()
        # MSFile
        # if self.isInterruptionRequested():
        #     return
        if fasta.exists():
            self.copy_file_to_folder(fasta, All_Cache)
        if rules.exists():
            self.copy_file_to_folder(rules, All_Cache)
        PAtO = []
        pool = QThreadPool.globalInstance()
        pool.setMaxThreadCount(self.containerThread)
        for group in P:
            PAt = []
            for path in group:
                path_obj = pathlib.Path(path) if not isinstance(path, pathlib.Path) else path
                if path_obj.exists() and path_obj.is_file():
                    PAt.append(path_obj.name)
                    self.copy_file_to_folder(path_obj, All_Cache)
                elif not path_obj.exists():
                    raise FileNotFoundError(f'{path_obj} does not exist')
                elif not path_obj.is_dir():
                    raise ValueError(f'{path_obj} is not a file path')
                # (path, save, All_Cache)
                task = GlycProcessingTask(path, save, All_Cache, self.update_signal)
                pool.start(task)
                # exit_code = container.wait()
                # logging.info(f"Container {container.id} finished with exit code {exit_code}")
                # if exit_code == 0:
                #     logging.info("Command executed successfully, proceeding with next steps or files...")
                #     container.stop()
                # else:
                #     logging.error("Command execution failed, check the container logs for details.")
            PAtO.append(PAt)
            pool.waitForDone()
        # GlycReSoft
        # Input pre-process
        absolute_path1 = []
        absolute_path2 = []
        for g_index,group in enumerate(PAtO):
            for name in group:
                pattern = re.compile(fr'^Cache-glycopepitides-in-samples_{name}-glycopeptides\.csv$')
                for file_path in save.iterdir():
                    if file_path.is_file() and pattern.match(str(file_path.name)):
                        if g_index == 0:
                            absolute_path1.append(file_path)
                        elif g_index == 1:
                            absolute_path2.append(file_path)
        # Change the result files to adjust the path for analysis
        result_save = pathlib.Path(self.save_path) / 'GlycReSoft_Analysis_Result'
        glyc_reload = GlyReSoft(absolute_path1, absolute_path2, fasta)
        glyc_reload.run()
        glyc_analysis = glyc_analysis_process(absolute_path1, absolute_path2, result_save, glyc_reload.total_output, self.FcValue)
        Ratio3, final = glyc_analysis.run()
        # 保存路径，拼接文件名
        file_path = f'{result_save}\\GlycReSoft_Tem_file_for_FinalAnalysis.pkl'
        with open(file_path,'wb') as f:
            pickle.dump(final,f)

        logging.info(f'Result:{Ratio3}')
        self.update_signal.emit(f"{result_save}\\Protein_Ratio_Difference.png")
        self.debugging_signal.emit('Finished')
        # return "Finish"

    def handle_metamorpheus_task(self):
        logging.info(f'Start handle Data from MetaMorpheus')
        # Start Processing data in Docker
        print('Start handle data with O-Pair')
        script_root_path = pathlib.Path(__file__).resolve().parent
        All_Cache = script_root_path / "All_Cache"
        # self.debugging_signal.emit("Start dealing with results from Glyc")
        # Start Processing data in Docker
        P = [self.file_raw1, self.file_raw2]
        # save = pathlib.Path(self.save_path) / 'GlycReSoft_Result'
        # fasta = pathlib.Path(resource_path('resources')) / 'database_Glyc_Tem.fasta'
        # rules = pathlib.Path(fasta).parent / "rules_file.txt"

        save = pathlib.Path(self.save_path) / 'O-Pair_Result'
        fasta = pathlib.Path(resource_path('resources')) / 'database_O-Pair_Tem.fasta'
        toml = pathlib.Path(fasta).parent / "Cache.toml"
        Setting = pathlib.Path(fasta).parent / 'MetaMorpheus'
        while not self.isInterruptionRequested():
            if len(self.MSF_fasta_output) == 0:
                self.debugging_signal.emit('Something unexpected has happened, please check the MS Files input. Task interrupted.')
                self.requestInterruption()
                break
            else:
                with open(fasta, 'w') as dst_file:
                    for sequence in self.MSF_fasta_output:
                        dst_file.write(sequence + "\n")
                break
        if fasta.exists():
            self.copy_file_to_folder(fasta, All_Cache)
        if toml.exists():
            self.copy_file_to_folder(toml, All_Cache)
        pool = QThreadPool.globalInstance()
        pool.setMaxThreadCount(self.containerThread)
        for group in P:
            for path in group:
                path_obj = pathlib.Path(path) if not isinstance(path, pathlib.Path) else path
                if path_obj.exists() and path_obj.is_file():
                    self.copy_file_to_folder(path_obj, All_Cache)
                elif not path_obj.exists():
                    raise FileNotFoundError(f'{path_obj} does not exist')
                elif not path_obj.is_dir():
                    raise ValueError(f'{path_obj} is not a file path')
                save_inner = save / path_obj.name
                task = MetaProcessingTask(path, save_inner, All_Cache, Setting, self.update_signal)
                pool.start(task)
            pool.waitForDone()
        # Design the path of psm and psmtsv files
        name_g = {}
        for index, group in enumerate(P):
            name = []
            for path in group:
                name.append(path)
            name_g[index] = name
        psm_g = {}
        peaks_g = {}
        for index, group in name_g.items():
            psm_g[index] = []
            peaks_g[index] = []
            for path in group:
                for folder in save.glob("*"):
                    if type(path) == pathlib.Path:
                        if folder.is_dir() and re.search(path.name, str(folder)):
                            all_psm = folder / "Task1GlycoSearchTask" / "oglyco.psmtsv"
                            psm_g[index].append(all_psm)
                            AllQuantifiedPeaks = folder / "Task1GlycoSearchTask" / "AllQuantifiedPeaks.tsv"
                            peaks_g[index].append(AllQuantifiedPeaks)
                    else:
                        raise TypeError(f'{path} is not a file path, something unexpected has happened')
        psm1 = psm_g[0]
        psm2 = psm_g[1]
        peaks1 = peaks_g[0]
        peaks2 = peaks_g[1]
        # Return to the normal process with Original Codes Part
        save_path = pathlib.Path(self.save_path) / "O-Pair_Result_Analysis"
        MetaReload = MetaMorpheus_Fasta_Reload(psm1, psm2, fasta, save_path)
        Me_R_O = MetaReload.run()
        MetaMorpheus = MetaMorpheusAnalysisProcess(psm1, psm2, peaks1, peaks2, save_path, Me_R_O, self.FcValue)
        Ratio3, final = MetaMorpheus.run()
        # 保存路径，拼接文件名
        file_path = f'{save_path}\\MetaMorpheus_Tem_file_for_FinalAnalysis.pkl'
        with open(file_path, 'wb') as f:
            pickle.dump(final, f)
        logging.info(f'Result:{Ratio3}')
        self.update_signal.emit(f"{save_path}\\Protein_Ratio_Difference.png")
        self.debugging_signal.emit('Finish')

    def handle_final_analysis(self):
        logging.info(f'Start handle Data from the whole Analysis Process')
        PCA_Glycosylation_Analysis.PCA_Glycosylation_Analysis(self.file_raw1, self.file_raw2, self.file_psm1, self.file_psm2, self.save_path,self.Annotation)

    def start_Statistic(self, Ratio, save_path):
        tem1 = []
        tem2 = []
        for key, value in Ratio.items():
            if key[0] == 1:
                tem1.append(value)
            elif key[0] == 2:
                tem2.append(value)
        Inner_Both_tem = []
        for value1 in tem1:
            for value2 in value1[0]:
                Inner_Both_tem.append(value2)
        for value1 in tem2:
            for value2 in value1[0]:
                Inner_Both_tem.append(value2)

        Inner_Both = []
        for item in Inner_Both_tem:
            if Inner_Both_tem.count(item) >= 6:
                Inner_Both.append(item)
        Inner_Both = list(set(Inner_Both))
        results_1 = []
        results_2 = []
        for key3, item in enumerate(tem1):
            for key4, value4 in enumerate(item[0]):
                for value3 in Inner_Both:
                    if value4 == value3:
                        r = item[1][key4]
                        ID = value4
                        r1 = (key3, ID, r)
                        results_1.append(r1)
        for key3, item in enumerate(tem2):
            for key4, value4 in enumerate(item[0]):
                for value3 in Inner_Both:
                    if value4 == value3:
                        r = item[1][key4]
                        ID = value4
                        r1 = (key3, ID, r)
                        results_2.append(r1)
        final_x = {}
        final_y = {}
        for item in Inner_Both:
            x = []
            for result in results_1:
                if item == result[1]:
                    x.append(result[2])
            final_x[item] = x
        for item in Inner_Both:
            y = []
            for result in results_2:
                if item == result[1]:
                    y.append(result[2])
            final_y[item] = y
        Sta = DataProvider(final_x, final_y, save_path, self.FcValue, self.PValue)
        log2fc, p_values, labels = Sta.generate_data()
        self.volcano_data_ready.emit(log2fc, p_values, labels, self.FcValue, self.PValue)

class MainApp(QMainWindow, Ui_MainWindow):
    debugging_signal = pyqtSignal(str)

    def __init__(self):
        super(MainApp, self).__init__()
        self.setupUi(self)
        # self.textBrowser = self.textBrowser
        # self.file_logger = Logger.Logger(level='debug')
        # print("Logger initialized successfully.")
        self.setFixedSize(self.size())
        self.config = load_config(CONFIG_FILE)
        self.graphicsView.setScene(QGraphicsScene())
        self.analysis_thread = None
        self.MSFragger_path = self.config.get('MSFragger_path', '')
        self.GlycReSoft_path = self.config.get('GlycReSoft_path', '')
        self.MetaMorpheus_path = self.config.get('MetaMorpheus_path', '')

        self.debugging_signal.connect(self.debugging_signal_show)

        # self.radio_group = QButtonGroup(self)
        # self.radio_group1 = QButtonGroup(self)
        self.output_signal = Signal()
        self.filePaths = {}
        self.pushButton.clicked.connect(self.fasta1_path)

        self.actionSetting.triggered.connect(self.Docker_Threads_Setting)
        self.pushButton_3.clicked.connect(self.get_path3)
        self.pushButton_4.clicked.connect(self.get_path4)
        self.pushButton_5.clicked.connect(self.select_directory)
        self.pushButton_6.clicked.connect(self.start_processing_data)
        self.pushButton_7.clicked.connect(self.exit_application)
        self.pushButton_8.clicked.connect(self.get_path5)
        self.pushButton_9.clicked.connect(self.get_path6)

        self.pushButton_14.clicked.connect(self.up_3)
        self.pushButton_15.clicked.connect(self.down_3)
        self.pushButton_16.clicked.connect(self.up_4)
        self.pushButton_17.clicked.connect(self.down_4)
        self.pushButton_18.clicked.connect(self.up_5)
        self.pushButton_19.clicked.connect(self.down_5)
        self.pushButton_20.clicked.connect(self.up_6)
        self.pushButton_21.clicked.connect(self.down_6)

        self.pushButton_27.clicked.connect(self.Delete_All_3)
        self.pushButton_29.clicked.connect(self.Delete_All_4)
        self.pushButton_31.clicked.connect(self.Delete_All_5)
        self.pushButton_33.clicked.connect(self.Delete_All_6)

        self.pushButton_26.clicked.connect(lambda: self.remove_selected_items(self.listWidget_3))
        self.pushButton_28.clicked.connect(lambda: self.remove_selected_items(self.listWidget_4))
        self.pushButton_30.clicked.connect(lambda: self.remove_selected_items(self.listWidget_5))
        self.pushButton_32.clicked.connect(lambda: self.remove_selected_items(self.listWidget_6))


        self.actionN_Glycosylation.setCheckable(True)
        self.actionO_Glycosylation.setCheckable(True)
        self.actionStep2.setCheckable(True)
        self.action_group = QtWidgets.QActionGroup(self)
        self.action_group.setExclusive(True)

        self.action_group.addAction(self.actionN_Glycosylation)
        self.action_group.addAction(self.actionO_Glycosylation)
        self.action_group.addAction(self.actionStep2)
        self.choice = 1

        # self.analysis_signal = MSFraggerRawRead
        # self.analysis_signal.text_update.connect(self.update_text_browser)
        #
        self.AnalysisSignal = MSFragger_raw_reading_signal
        self.AnalysisSignal.text_update.connect(self.update_text_browser)

        self.outputSignal = MSFragger_reload_signal
        self.outputSignal.text_update.connect(self.update_text_browser)

        # self.debug_signal.connect(self.update_text_browser)
        self.action_group.triggered.connect(self.choose)

        self.actionAbout.triggered.connect(self.open_about)
        self.actionOpen_the_Save_Folder.triggered.connect(self.open_save_folder)
        self.actionOutput_the_Log.triggered.connect(self.output_log)
        self.actionPath_Set.triggered.connect(self.set_MSFragger_path)
        self.actionStart_3.triggered.connect(self.open_MSFragger)
        self.actionPath_Set_2.triggered.connect(self.set_GlycReSoft_path)
        self.actionStart.triggered.connect(self.open_GlycReSoft)
        self.actionPath_Set_3.triggered.connect(self.set_MetaMorpheus_path)
        self.actionStart_2.triggered.connect(self.open_MetaMorpheus)

    def get_current_log_path(self):
        return self.file_logger.get_log_file()

    @pyqtSlot(str)
    def update_text_browser(self, text):
        cursor = self.textBrowser.textCursor()
        cursor.movePosition(QTextCursor.End)
        self.textBrowser.append(text)
        self.textBrowser.setTextCursor(cursor)
        self.textBrowser.ensureCursorVisible()

    @pyqtSlot(list, list, list, float, float)
    def display_volcano_plot(self, log2fc, p_values, labels, fcvalue, pvalue):
        self.volcano_window = ImageWindow(log2fc, p_values, labels, fcvalue, pvalue)
        self.volcano_window.show()

    @pyqtSlot(str)
    def debugging_signal_show(self, message):
        try:
            if "interrupt" in message:
                QMessageBox.information(self, 'Message', message, QMessageBox.Ok)
            self.analysis_thread.quit()
            # self.analysis_thread.wait()
        except Exception as e:
            print(e)


    def stop_thread(self):
        if self.analysis_thread and self.analysis_thread.isRunning():
            self.analysis_thread.quit()
            # self.analysis_thread.wait()
            self.update_text_browser('Stopping thread...')
            time.sleep(2)
            self.analysis_thread.stop()

    def analysis_finished(self):
        self.update_text_browser("Analysis Thread Finished.")
        # self.radioButton.setEnabled(True)
        # self.radioButton_2.setEnabled(True)
        # self.radioButton_3.setEnabled(True)

    # def showComboBox(self):
    #     items = ["N-Glycopeptide", "O-Glycopeptide"]
    #     item, ok = QInputDialog.getItem(self, "Choose", "Please Choose the Running Mode:", items, 0, False)
    #     if ok and item:
    #         return item

    def start_processing_data(self):
        self.pushButton_6.setEnabled(False)
        if self.actionO_Glycosylation.isChecked() == True:
            value_2 = 'O-Glycopeptide'
            value_1, ok1 = QInputDialog.getDouble(self, "Fold Change", "\nInput Decimal:", 1.50, 0.00, 10000.00, 2)
            value_3, ok3 = QInputDialog.getDouble(self, "P-Value", "\nInput Decimal:", 0.05, 0.00, 10000.00, 2)
            value_ensure = QMessageBox.question(self,'Notice', " Please ensure the setting of MetaMorpheus is correct", QMessageBox.Yes | QMessageBox.No)

            if value_ensure == QMessageBox.Yes:
                pass
            else:
                self.pushButton_6.setEnabled(True)
                return
            if ok1 and ok3:
                value1 = value_1  # FC阈值
                value2 = value_3
                # self.pushButton_6.setEnabled(False)
                try:
                    self.filePaths['RAW1'] = [self.listWidget_3.item(i).text() for i in range(self.listWidget_3.count())]
                    self.filePaths['RAW2'] = [self.listWidget_5.item(i).text() for i in range(self.listWidget_5.count())]
                    self.filePaths['PSM1'] = [self.listWidget_4.item(i).text() for i in range(self.listWidget_4.count())]
                    self.filePaths['PSM2'] = [self.listWidget_6.item(i).text() for i in range(self.listWidget_6.count())]
                except Exception as e:
                    QMessageBox.critical(self, "Warning", f"An error occurred: {e}", QMessageBox.Ok)
                    self.pushButton_6.setEnabled(True)
                file_raw1 = self.filePaths['RAW1']
                file_raw2 = self.filePaths['RAW2']
                file_peak1 = self.filePaths['PSM1']
                file_peak2 = self.filePaths['PSM2']
                save_path = self.filePaths['SavePath']
                for file_path in [file_raw1, file_raw2, file_peak1, file_peak2, save_path]:
                    if isinstance(file_path, list):
                        for path in file_path:
                            if not os.path.exists(path):
                                self.update_text_browser(f"Cannot find: {path}")
                                self.pushButton_6.setEnabled(True)
                                return
                    else:
                        if not os.path.exists(file_path):
                            self.update_text_browser(f"Cannot find: {file_path}")
                            self.pushButton_6.setEnabled(True)
                            return
                try:
                    # file_psm1, file_psm2, save_path, output_signal, analysis_signal, raw_reading_signal, reload_signal,
                    # file_raw1, file_raw2, fasta1_filepath, choice, value, value_2, FcValue, PValue
                    self.analysis_thread = AnalysisThread(
                        file_raw1, file_raw2, file_peak1, file_peak2, save_path,
                        self.output_signal, self.analysis_signal,
                        self.outputSignal, fasta1_filepath=self.fasta1_filepath,
                        choice=2, value_2=value_2, FcValue=value1, PValue=value2,ContainerThread = self.Container_Thread_int,Annotation=None
                    )
                    self.analysis_thread.update_signal.connect(self.update_plot)
                    self.analysis_thread.volcano_data_ready.connect(self.display_volcano_plot)
                    self.analysis_thread.finished.connect(lambda: self.pushButton_6.setEnabled(True))
                    self.analysis_thread.start()
                except Exception as e:
                    QMessageBox.critical(self, "Warning", f"An error occurred: {e}", QMessageBox.Ok)
                    self.pushButton_6.setEnabled(True)
                    # self.radioButton.setEnabled(True)
                    # self.pushButton_6.setEnabled(True)
                finally:
                    # self.radioButton.setEnabled(True)
                    self.listWidget_4.setEnabled(True)
                    self.listWidget_6.setEnabled(True)
                    self.pushButton_8.setEnabled(True)
                    self.pushButton_9.setEnabled(True)
                    self.pushButton_16.setEnabled(True)
                    self.pushButton_17.setEnabled(True)
                    self.pushButton_20.setEnabled(True)
                    self.pushButton_21.setEnabled(True)
                    self.pushButton_28.setEnabled(True)
                    self.pushButton_29.setEnabled(True)
                    self.pushButton_32.setEnabled(True)
                    self.pushButton_33.setEnabled(True)
            else:
                QMessageBox.warning(self, "Warning", "Cancel Input")
                self.pushButton_6.setEnabled(True)
        elif self.actionN_Glycosylation.isChecked() == True:
            value_2 = 'N-Glycopeptide'
            value_1, ok1 = QInputDialog.getDouble(self, "Fold Change", "\nInput Decimal:", 1.50, 0.00, 10000.00, 2)
            value_3, ok3 = QInputDialog.getDouble(self, "P-Value", "\nInput Decimal:", 0.05, 0.00, 10000.00, 2)
            if ok1 and ok3:
                value1 = value_1  # FC阈值
                value2 = value_3  # P阈值
                self.filePaths['PSM1'] = [self.listWidget_3.item(i).text() for i in range(self.listWidget_3.count())]
                self.filePaths['PSM2'] = [self.listWidget_5.item(i).text() for i in range(self.listWidget_5.count())]
                self.filePaths['Peak1'] = [self.listWidget_4.item(i).text() for i in range(self.listWidget_4.count())]
                self.filePaths['Peak2'] = [self.listWidget_6.item(i).text() for i in range(self.listWidget_6.count())]

                if isinstance(value1, (float, int)):
                    try:
                        file_raw1 = self.filePaths['PSM1']
                        file_raw2 = self.filePaths['PSM2']
                        Total = self.filePaths['Peak1']
                        Glycosylation = self.filePaths['Peak2']
                        save_path = self.filePaths['SavePath']
                        fasta1_filepath = self.fasta1_filepath

                        for file_path in [file_raw1, file_raw2, Total, Glycosylation, save_path, fasta1_filepath]:
                            if isinstance(file_path, list):
                                for path in file_path:
                                    if not os.path.exists(path):
                                        self.update_text_browser(f"Cannot find: {path}")
                                        self.pushButton_6.setEnabled(True)
                                        return
                            else:
                                if not os.path.exists(file_path):
                                    self.update_text_browser(f"Cannot find: {file_path}")
                                    self.pushButton_6.setEnabled(True)
                                    return

                        try: # file_raw1, file_raw2, peak1, peak2, save_path, output_signal, analysis_signal, reload_signal, fasta1_filepath, choice, value, value_2, FcValue, PValue
                            self.analysis_thread = AnalysisThread(
                                file_raw1, file_raw2, Total, Glycosylation, save_path,
                                self.output_signal, self.analysis_signal,
                                self.outputSignal, fasta1_filepath,
                                self.choice, value_2, value1, value2, ContainerThread = self.Container_Thread_int,Annotation=None # P阈值
                            )   # N/O Bottom  N/O MSF   FC     P
                            self.analysis_thread.debugging_signal.connect(self.debugging_signal_show)
                            self.analysis_thread.update_signal.connect(self.update_plot)
                            self.analysis_thread.volcano_data_ready.connect(self.display_volcano_plot)
                            self.analysis_thread.finished.connect(lambda: self.pushButton_6.setEnabled(True))
                            self.analysis_thread.start()
                        except Exception as er:
                            print(er)
                            self.stop_thread()
                            self.pushButton_6.setEnabled(True)
                        finally:
                            # self.radioButton.setEnabled(True)
                            self.listWidget_4.setEnabled(True)
                            self.listWidget_6.setEnabled(True)
                            self.pushButton_8.setEnabled(True)
                            self.pushButton_9.setEnabled(True)
                            self.pushButton_16.setEnabled(True)
                            self.pushButton_17.setEnabled(True)
                            self.pushButton_20.setEnabled(True)
                            self.pushButton_21.setEnabled(True)
                            self.pushButton_28.setEnabled(True)
                            self.pushButton_29.setEnabled(True)
                            self.pushButton_32.setEnabled(True)
                            self.pushButton_33.setEnabled(True)

                        self.analysis_thread.finished.connect(self.analysis_finished)
                    except Exception as e:
                        QMessageBox.critical(self, "Warning", f"An error occurred: {e} during N-glycosylation Analysis", QMessageBox.Ok)
                        self.stop_thread()
                        self.pushButton_6.setEnabled(True)
                        # self.pushButton_6.setEnabled(True)
                    # finally:
                    #     self.listWidget_4.setEnabled(True)
                    #     self.listWidget_6.setEnabled(True)
                    #     self.pushButton_8.setEnabled(True)
                    #     self.pushButton_9.setEnabled(True)
                    #     self.pushButton_16.setEnabled(True)
                    #     self.pushButton_17.setEnabled(True)
                    #     self.pushButton_20.setEnabled(True)
                    #     self.pushButton_21.setEnabled(True)
                    #     self.pushButton_28.setEnabled(True)
                    #     self.pushButton_29.setEnabled(True)
                    #     self.pushButton_32.setEnabled(True)
                    #     self.pushButton_33.setEnabled(True)
                else:
                    QMessageBox.critical(self, "Warning", "An error occurred: Choose the right value", QMessageBox.Ok)
                    # self.stop_thread()
                    self.pushButton_6.setEnabled(True)
            else:
                QMessageBox.warning(self, "Warning", "Cancel Input")
                self.pushButton_6.setEnabled(True)
        else:
            # t = tk.Tk()
            # t.withdraw()
            judge = QMessageBox.question(
                self,
                'Notice',
                'Requires input of raw mass spectrometry files',
                QMessageBox.Yes | QMessageBox.No  # 更常见的组合
            )
            Annotation = None
            if judge == QMessageBox.Yes:
                filegroup1, _ = QFileDialog.getOpenFileNames(self,'Open file for Experimental group','C:\\','MS files (*.raw *.mzml)')
                filegroup2, _ = QFileDialog.getOpenFileNames(self, 'Open file for Control group', 'C:\\','MS files (*.raw *.mzml)')
                if not filegroup2 or not filegroup1:
                    QMessageBox.critical(self, "Warning",
                                         f"An error occurred: Fatel error during setting MS filepaths, check the parameters",
                                         QMessageBox.Ok)
                    # self.stop_thread()
                    self.pushButton_6.setEnabled(True)
                else:
                    Annotation = [filegroup1,filegroup2]

            else:
                QMessageBox.warning(self, "Warning", "Cancel Input")
                self.pushButton_6.setEnabled(True)
                return

            self.filePaths['Excel1'] = [self.listWidget_3.item(i).text() for i in range(self.listWidget_3.count())]
            self.filePaths['Excel2'] = [self.listWidget_5.item(i).text() for i in range(self.listWidget_5.count())]
            self.filePaths['Total'] = [self.listWidget_4.item(i).text() for i in range(self.listWidget_4.count())]
            self.filePaths['Glycosylation'] = [self.listWidget_6.item(i).text() for i in range(self.listWidget_6.count())]
            file_raw1 = None
            file_raw2 = None
            Total = None
            Glycosylation = None
            save_path = None
            try:
                file_raw1 = self.filePaths['Excel1']
                file_raw2 = self.filePaths['Excel2']
                Total = self.filePaths['Total']
                Glycosylation = self.filePaths['Glycosylation']
                save_path = self.filePaths['SavePath']
                for file_path in [file_raw1, file_raw2, Total, Glycosylation, save_path]:
                    if isinstance(file_path, list):
                        for path in file_path:
                            if not os.path.exists(path):
                                self.update_text_browser(f"Cannot find: {path}")
                                self.pushButton_6.setEnabled(True)
                                return
                    else:
                        if not os.path.exists(file_path):
                            self.update_text_browser(f"Cannot find: {file_path}")
                            self.pushButton_6.setEnabled(True)
                            return
            except Exception as e:
                print(e)
                self.stop_thread()
                self.pushButton_6.setEnabled(True)
            # finally:
            #     self.pushButton_6.setEnabled(True)
            try:  # file_raw1, file_raw2, file_psm1, file_psm2, save_path, output_signal, analysis_signal, reload_signal, fasta1_filepath, choice, value, value_2, FcValue, PValue
                self.analysis_thread = AnalysisThread(
                    file_raw1=file_raw1, file_raw2=file_raw2, file_psm1=Total, file_psm2=Glycosylation, save_path=save_path,
                    output_signal=None, analysis_signal=self.analysis_signal,
                    reload_signal=None, fasta1_filepath=None,
                    choice=None, value_2=None, FcValue=None, PValue=None,ContainerThread=1,Annotation=Annotation  # P阈值
                )  # N/O Bottom  N/O MSF   FC     P
                self.analysis_thread.debugging_signal.connect(self.debugging_signal_show)
                self.analysis_thread.update_signal.connect(self.update_plot)
                # self.analysis_thread.volcano_data_ready.connect(self.display_volcano_plot)
                self.analysis_thread.finished.connect(lambda: self.pushButton_6.setEnabled(True))
                self.analysis_thread.start()
            except Exception as er:
                print(er)
                self.stop_thread()
                self.pushButton_6.setEnabled(True)
            # finally:
            #     self.stop_thread()
            #     self.pushButton_6.setEnabled(True)

    @pyqtSlot(str)
    def update_plot(self, filepath):
        logging.info(f"Received update signal with filepath: {filepath}")
        if os.path.exists(filepath):
            self.graphicsView.scene().clear()
            pixmap = QtGui.QPixmap(filepath)
            self.pixmap_item = QGraphicsPixmapItem(pixmap)
            self.graphicsView.scene().addItem(self.pixmap_item)
            logging.info(f"Image updated in QGraphicsView with filepath: {filepath}")
            self.update_text_browser("update of Graphics")

            # 使图像适应视图窗口大小
            self.graphicsView.fitInView(self.pixmap_item, Qt.KeepAspectRatio)
            # 将视图中心移动到图像中心
            self.graphicsView.centerOn(self.pixmap_item)

            # 添加鼠标滚轮事件以支持缩放
            self.graphicsView.viewport().installEventFilter(self)

        else:
            self.update_text_browser("Without update of Graphics")
            logging.error(f"File not found: {filepath}")

    def eventFilter(self, source, event):
        if event.type() == QtCore.QEvent.Wheel and source is self.graphicsView.viewport():
            if event.modifiers() == Qt.ControlModifier:
                if event.angleDelta().y() > 0:
                    self.graphicsView.scale(1.2, 1.2)  # 放大
                else:
                    self.graphicsView.scale(0.8, 0.8)  # 缩小
                return True
        return super(MainApp, self).eventFilter(source, event)
    @classmethod
    def move_up(self, path_list_widget):
        current_row = path_list_widget.currentRow()
        if current_row > 0:
            item = path_list_widget.takeItem(current_row)
            path_list_widget.insertItem(current_row - 1, item)
            path_list_widget.setCurrentRow(current_row - 1)

    @classmethod
    def move_down(self, path_list_widget):
        current_row = path_list_widget.currentRow()
        if current_row < path_list_widget.count() - 1:
            item = path_list_widget.takeItem(current_row)
            path_list_widget.insertItem(current_row + 1, item)
            path_list_widget.setCurrentRow(current_row + 1)

    def up_3(self):
        self.move_up(self.listWidget_3)

    def down_3(self):
        self.move_down(self.listWidget_3)

    def up_4(self):
        self.move_up(self.listWidget_4)

    def down_4(self):
        self.move_down(self.listWidget_4)

    def up_5(self):
        self.move_up(self.listWidget_5)

    def down_5(self):
        self.move_down(self.listWidget_5)

    def up_6(self):
        self.move_up(self.listWidget_6)

    def down_6(self):
        self.move_down(self.listWidget_6)

    def Delete_All_3(self):
        self.listWidget_3.clear()

    def Delete_All_4(self):
        self.listWidget_4.clear()

    def Delete_All_5(self):
        self.listWidget_5.clear()

    def Delete_All_6(self):
        self.listWidget_6.clear()

    def remove_selected_items(self, list_widget):
        selected_items = list_widget.selectedItems()
        for item in selected_items:
            list_widget.takeItem(list_widget.row(item))

    def fasta1_path(self):
        path14, _ = QFileDialog.getOpenFileName(self, "Choose File", "", "FASTA Files (*.fasta *.fa *.fas)")
        self.fasta1_filepath = path14
        self.lineEdit.setText(path14)

    def get_path3(self):
        if self.actionStep2.isChecked() == True:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Byonic Output Experiment Group Files", "", "xlsx Files (*.xlsx)")
            if selected_files:
                for file in selected_files:
                    self.listWidget_3.addItem(file)
        # elif self.radioButton_4.isChecked() == True:
        #     selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Group I Mass Spectrometry Files", "","mzML Files (*.mzml)")
        #     if selected_files:
        #         self.listWidget_3.addItems(selected_files)
        elif self.actionO_Glycosylation.isChecked() == True:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Group I Mass Spectrometry Files", "", "Raw Files (*.raw)")
            if selected_files:
                for file in selected_files:
                    self.listWidget_3.addItem(file)
        else:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Group I Mass Spectrometry Files", "", "MZML Files (*.mzml)")
            if selected_files:
                for file in selected_files:
                    self.listWidget_3.addItem(file)

    def get_path4(self):
        if self.actionStep2.isChecked() == True:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Byonic Output Control Group Files", "", "xlsx Files (*.xlsx)")
            if selected_files:
                for file in selected_files:
                    self.listWidget_5.addItem(file)
        # elif self.radioButton_4.isChecked() == True:
        #     selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Group II Mass Spectrometry Files", "","mzML Files (*.mzml)")
        #     if selected_files:
        #         self.listWidget_5.addItems(selected_files)
        elif self.actionO_Glycosylation.isChecked() == True:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Group II Mass Spectrometry Files", "", "Raw Files (*.raw)")
            if selected_files:
                for file in selected_files:
                    self.listWidget_5.addItem(file)
        else:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Group II Mass Spectrometry Files", "", "mzML Files (*.mzml)")
            if selected_files:
                for file in selected_files:
                    self.listWidget_5.addItem(file)

    def get_path5(self):
        if self.actionStep2.isChecked() == True:
            directory = QFileDialog.getExistingDirectory(self, "Choose Folder for Total Analysis output folder", "./", QFileDialog.ShowDirsOnly)
            self.filePaths["T"] = directory.replace('/', '\\')
            self.listWidget_4.addItem(directory)
        else:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Experiment Group result from MSFragger", "", 'TSV Files (*.tsv)')
            if selected_files:
                for file in selected_files:
                    self.listWidget_4.addItem(file)

    def get_path6(self):
        if self.actionStep2.isChecked() == True:
            directory = QFileDialog.getExistingDirectory(self, "Choose Folder for N/O-Glycosylation Analysis output folder", "./", QFileDialog.ShowDirsOnly)
            self.filePaths["A"] = directory.replace('/', '\\')
            self.listWidget_6.addItem(directory)
        else:
            selected_files, _ = QFileDialog.getOpenFileNames(self, "Choose Control Group result from MSFragger", "", 'TSV Files (*.tsv)')
            if selected_files:
                for file in selected_files:
                    self.listWidget_6.addItem(file)

    def select_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Folder", "./", QFileDialog.ShowDirsOnly)
        self.lineEdit_5.setText(directory)
        self.filePaths["SavePath"] = directory.replace('/', '\\')
        self.update_text_browser(f'Setting the Save Path: {self.filePaths["SavePath"]}')

    def exit_application(self):
        self.stop_thread()
        sys.exit()


    def open_about(self):
        subprocess.Popen(['notepad', resource_path('resources/About')])

    def Docker_Threads_Setting(self):
        value, ok = QInputDialog.getInt(self, "Threads Number", "\nInput Decimal:", 1, 1, 3, 1)
        if ok:
            self.Container_Thread_int = value
        else:
            QMessageBox.information(self, "Thread Notice", "The program will use the default thread number to activate task pool in Docker")
    def open_save_folder(self):
        try:
            directory = self.filePaths["SavePath"]
            if not directory:
                QMessageBox.critical(self, "Warning", f"SavePath folder path is not set", QMessageBox.Ok)
                return
            if QtCore.QSysInfo.productType() == "windows":
                subprocess.Popen(['explorer', directory], shell=True)
            elif QtCore.QSysInfo.productType() == "osx":
                subprocess.Popen(["open", directory])
            elif QtCore.QSysInfo.productType() == "linux":
                subprocess.Popen(["xdg-open", directory])
        except Exception as e:
            QMessageBox.critical(self, 'Error', f"Error opening file explorer:{e}", QMessageBox.Ok)

    def set_MSFragger_path(self):
        path10, _ = QFileDialog.getOpenFileName(self, "Choose File", "", "EXE Files (*.exe)")
        if path10:
            self.MSFragger_path = path10
            self.update_text_browser(f'Setting MSFragger Path:{path10}')
            self.config['MSFragger_path'] = self.MSFragger_path
            save_config(self.config)

    def open_MSFragger(self):
        try:
            exe_path = self.MSFragger_path
            if not exe_path:
                QMessageBox.critical(self, "Warning", "MetaMorpheus path is not set", QMessageBox.Ok)
                return
            if QtCore.QSysInfo.productType() == "windows":
                subprocess.Popen([exe_path], shell=True)
            else:
                QMessageBox.critical(self, "Warning", "Unsupported OS", QMessageBox.Ok)
        except Exception as e:
            QMessageBox.critical(self, "Warning", f"Error opening MetaMorpheus: {e}", QMessageBox.Ok)

    def set_GlycReSoft_path(self):
        path11, _ = QFileDialog.getOpenFileName(self, "Choose File", "", "EXE Files (*.exe)")
        if path11:
            self.GlycReSoft_path = path11
            self.update_text_browser(f'Setting GlycReSoft Path:{path11}')
            self.config['GlycReSoft_path'] = self.GlycReSoft_path
            save_config(self.config)

    def open_GlycReSoft(self):
        try:
            exe_path = self.GlycReSoft_path
            if not exe_path:
                QMessageBox.critical(self, "Warning", "MetaMorpheus path is not set", QMessageBox.Ok)
                return
            if QtCore.QSysInfo.productType() == "windows":
                subprocess.Popen([exe_path], shell=True)
            else:
                QMessageBox.critical(self, "Warning", "Unsupported OS", QMessageBox.Ok)
        except Exception as e:
            QMessageBox.critical(self, "Warning", f"Error opening MetaMorpheus: {e}", QMessageBox.Ok)

    def set_MetaMorpheus_path(self):
        path12, _ = QFileDialog.getOpenFileName(self, "Choose File", "", "EXE Files (*.exe)")
        if path12:
            self.MetaMorpheus_path = path12
            self.update_text_browser(f'Setting MetaMorpheus Path:{path12}')
            self.config['MetaMorpheus_path'] = self.MetaMorpheus_path
            save_config(self.config)

    def open_MetaMorpheus(self):
        try:
            exe_path = self.MetaMorpheus_path
            if not exe_path:
                QMessageBox.critical(self, "Warning", "MetaMorpheus path is not set", QMessageBox.Ok)
                return
            if QtCore.QSysInfo.productType() == "windows":
                subprocess.Popen([exe_path], shell=True)
            else:
                QMessageBox.critical(self, "Warning", "Unsupported OS", QMessageBox.Ok)
        except Exception as e:
            QMessageBox.critical(self, "Warning", f"Error opening MetaMorpheus: {e}", QMessageBox.Ok)

    def output_log(self):
        log_file = self.get_current_log_path()
        day = time.strftime('%Y-%m-%d_%H-%M-%S', time.localtime())
        try:
            dest_folder = QFileDialog.getExistingDirectory(self, "Choose Folder", "./", QFileDialog.ShowDirsOnly)
            if not dest_folder:
                return

            dest_file = os.path.join(dest_folder, f"{day}_exported.log")
            with open(log_file, 'r', encoding='utf-8') as src, open(dest_file, 'w', encoding='utf-8') as dst:
                dst.writelines(src.readlines())

            QMessageBox.information(self, "Success", "Log File Exported Successfully!", QMessageBox.Ok)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error Exporting Log File: {e}", QMessageBox.Ok)

    @pyqtSlot()
    def choose(self):
        try:
            self.output_signal.text_update.disconnect()
            self.analysis_signal.text_update.disconnect()
            # self.AnalysisSignal.text_update.disconnect()
        except TypeError:
            pass
        # if self.btn1.isChecked():
        #     self.choice = 1
        #     self.label_4.setText('Choose Group I Mass Spectrometry Files')
        #     self.label_5.setText('Choose Group II Mass Spectrometry Files')
        #     self.output_signal = MSFragger_reload_signal
        #     # self.analysis_signal = MSFraggerRawRead
        #     self.analysis_signal = MSFragger_raw_reading_signal
        #     self.lineEdit.setEnabled(True)
        #     self.pushButton.setEnabled(True)
        #     self.listWidget_4.setEnabled(False)
        #     self.listWidget_6.setEnabled(False)
        #     self.pushButton_8.setEnabled(False)
        #     self.pushButton_9.setEnabled(False)
        #     self.pushButton_16.setEnabled(False)
        #     self.pushButton_17.setEnabled(False)
        #     self.pushButton_20.setEnabled(False)
        #     self.pushButton_21.setEnabled(False)
        #     self.pushButton_28.setEnabled(False)
        #     self.pushButton_29.setEnabled(False)
        #     self.pushButton_32.setEnabled(False)
        #     self.pushButton_33.setEnabled(False)
        if self.actionN_Glycosylation.isChecked():
            self.choice = 1

            self.pushButton_3.setEnabled(True)
            self.pushButton_26.setEnabled(True)
            self.pushButton_27.setEnabled(True)
            self.pushButton_14.setEnabled(True)
            self.pushButton_15.setEnabled(True)

            self.pushButton_4.setEnabled(True)
            self.pushButton_18.setEnabled(True)
            self.pushButton_19.setEnabled(True)
            self.pushButton_30.setEnabled(True)
            self.pushButton_31.setEnabled(True)

            self.pushButton_6.setEnabled(True)

            self.lineEdit.setEnabled(True)
            self.pushButton.setEnabled(True)
            self.label_4.setText('Experimental Group MS Files Path')
            self.label_5.setText('Control Group MS Files Path')
            self.label_8.setText('MSFragger Result of Experimental Group')
            self.label_9.setText('MSFragger Result of Control Group')
            self.output_signal = GlyreloadOutput
            self.analysis_signal = GlyAnalysisoutput
            self.listWidget_4.setEnabled(True)
            self.listWidget_6.setEnabled(True)
            self.pushButton_8.setEnabled(True)
            self.pushButton_9.setEnabled(True)
            self.pushButton_16.setEnabled(True)
            self.pushButton_17.setEnabled(True)
            self.pushButton_20.setEnabled(True)
            self.pushButton_21.setEnabled(True)
            self.pushButton_28.setEnabled(True)
            self.pushButton_29.setEnabled(True)
            self.pushButton_32.setEnabled(True)
            self.pushButton_33.setEnabled(True)
        elif self.actionO_Glycosylation.isChecked():
            self.choice = 2
            self.pushButton_3.setEnabled(True)
            self.pushButton_26.setEnabled(True)
            self.pushButton_27.setEnabled(True)
            self.pushButton_14.setEnabled(True)
            self.pushButton_15.setEnabled(True)

            self.pushButton_4.setEnabled(True)
            self.pushButton_18.setEnabled(True)
            self.pushButton_19.setEnabled(True)
            self.pushButton_30.setEnabled(True)
            self.pushButton_31.setEnabled(True)

            self.pushButton_6.setEnabled(True)

            self.lineEdit.setEnabled(True)
            self.pushButton.setEnabled(True)
            self.label_4.setText('Experimental Group MS Files Path')
            self.label_5.setText('Control Group MS Files Path')
            self.label_8.setText('MSFragger Result of Experimental Group')
            self.label_9.setText('MSFragger Result of Control Group')
            self.output_signal = MetaReloadOutput
            self.analysis_signal = MetaMorpheusAnalysisoutput
            self.listWidget_4.setEnabled(True)
            self.listWidget_6.setEnabled(True)
            self.pushButton_8.setEnabled(True)
            self.pushButton_9.setEnabled(True)
            self.pushButton_16.setEnabled(True)
            self.pushButton_17.setEnabled(True)
            self.pushButton_20.setEnabled(True)
            self.pushButton_21.setEnabled(True)
            self.pushButton_28.setEnabled(True)
            self.pushButton_29.setEnabled(True)
            self.pushButton_32.setEnabled(True)
            self.pushButton_33.setEnabled(True)
        elif self.actionStep2.isChecked():
            self.choice = 3
            self.pushButton_3.setEnabled(True)
            self.pushButton_26.setEnabled(True)
            self.pushButton_27.setEnabled(True)
            self.pushButton_14.setEnabled(True)
            self.pushButton_15.setEnabled(True)

            self.pushButton_4.setEnabled(True)
            self.pushButton_18.setEnabled(True)
            self.pushButton_19.setEnabled(True)
            self.pushButton_30.setEnabled(True)
            self.pushButton_31.setEnabled(True)

            self.pushButton_6.setEnabled(True)

            self.label_4.setText('Byonic Analysis Result Experimental Group')
            self.label_5.setText('Byonic Analysis Result Control Group')
            self.label_8.setText('Total Analysis Output Folder')
            self.label_9.setText('N/O-Glycosylation Analysis Output Folder')
            self.output_signal = final_signal
            self.analysis_signal = final_signal
            self.listWidget_4.setEnabled(True)
            self.listWidget_6.setEnabled(True)
            self.pushButton_8.setEnabled(True)
            self.pushButton_9.setEnabled(True)
            self.pushButton_16.setEnabled(True)
            self.pushButton_17.setEnabled(True)
            self.pushButton_20.setEnabled(True)
            self.pushButton_21.setEnabled(True)
            self.pushButton_28.setEnabled(True)
            self.pushButton_29.setEnabled(True)
            self.pushButton_32.setEnabled(True)
            self.pushButton_33.setEnabled(True)
            self.lineEdit.setEnabled(False)
            self.pushButton.setEnabled(False)
        self.output_signal.text_update.connect(self.update_text_browser)

if __name__ == "__main__":
    try:
        QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
        QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)
        app = QApplication(sys.argv)
        win = MainApp()
        win.show()
        logging.debug("Starting application")
        app.exec_()
    except Exception as e:
        traceback.print_exc()
        print(f"Error opening file explorer:{e}")
    # finally:
    #