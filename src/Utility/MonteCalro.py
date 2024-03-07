import numpy as np
import time
class MonteBase:
    def __init__(self):
        self.start_time=0
        self.end_time=0
    def execution_time(self):
        return self.end_time-self.start_time
    def setSize(self,num):
        self.n=num
    def random(self,dd):
        #辞書型データに基づき、乱数を発生し戻す
        m=dd['mean']
        cov=dd['cov']
        s=m*cov
        generator = np.random.default_rng()
        if dd['dist']=='normal':
            ran = generator.normal(loc=m, scale=s, size=self.n)
        return ran
    def start(self):
        # 開始時間を記録
        self.start_time = time.time()
    def end(self):
        # 終了時間を記録
        self.end_time = time.time()
    def result(self):
        res={}
        res['n']=self.n
        res['Pf']=self.Pf
        res['Time']=self.execution_time()
        return res
    def monteG(self,data):
        #仮想関数
        pf=0
        print('class monteGが定義されていません')
        if type(data) is not dict:
            print('dataか辞書型ではありません')
        return pf
    def Calc(self,data):
        self.start()
        self.Pf=self.monteG(data)
        self.end()
    def setTitle(self,title,var):
        self.title=title
        self.var=var
    def getTitle(self):
        return self.title,self.var 