# Signal
from PyQt5.QtCore import QObject, pyqtSignal

class AnalysisSignal(QObject):
    text_update = pyqtSignal(str)
    # debugging_update = pyqtSignal(str)
    def write(self, text):
        self.text_update.emit(str(text))

    # def update(self, text):
    #     self.debugging_update.emit(str(text))