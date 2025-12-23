# Logger.py
# UTF-8

import time
import os
import logging
from logging import handlers

formater = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s'

def log_path_check():
    """Get the tem log file path"""
    log_path = os.path.join(os.getcwd(), "logs")
    if not os.path.exists(log_path):
        os.makedirs(log_path)
    timestamp = time.strftime('%Y-%m-%d_%H-%M-%S', time.localtime())
    log_name = f"temporary_{timestamp}.log"
    return os.path.join(log_path, log_name)


class Logger(object):
    def __init__(self, filename=log_path_check(), level='info', when='D', backCount=3, fmt=formater):
        self.level_relations = {
            'debug': logging.DEBUG,
            'info': logging.INFO,
            'warning': logging.WARNING,
            'error': logging.ERROR,
            'crit': logging.CRITICAL
        }
        self.logger = logging.getLogger(__name__)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()  # 清除已有的 Handler 防止重复日志
        format_str = logging.Formatter(fmt)
        self.logger.setLevel(self.level_relations.get(level))

        # 屏幕日志输出
        sh = logging.StreamHandler()
        sh.setFormatter(format_str)
        self.logger.addHandler(sh)

        # 文件日志输出（临时缓存）
        th = handlers.TimedRotatingFileHandler(filename=filename, when=when, backupCount=backCount, encoding='utf-8')
        th.setFormatter(format_str)
        self.logger.addHandler(th)

        self.filename = filename  # 暴露临时日志文件名

    def get_log_file(self):
        """获取当前日志文件路径"""
        return self.filename