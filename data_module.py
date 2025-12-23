# data_module.py

import T_test_Vocano
class DataProvider:
    def __init__(self, final_x, final_y, save_path, FcValue, PValue):
        self.final_x = final_x
        self.final_y = final_y
        self.save_path = save_path
        self.FcValue = FcValue
        self.PValue = PValue
        self.sta = T_test_Vocano.Statistic(self.final_x,self.final_y,self.save_path,self.FcValue,self.PValue)

    def generate_data(self):
        # 生成数据
        return self.sta.Volcano()

