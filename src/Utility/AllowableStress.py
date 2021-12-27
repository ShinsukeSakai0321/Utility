import statistics
import numpy as np
from scipy.stats import norm 
from scipy.stats import nct # 非心t分布
class PAStress:
    """
    目的:
        許容応力の確率論的決定支援(正規分布)
    例題:
        pas=PAStress()
        data=[408.4, 374.0, 395.9, 405.8,412.3]
        pas.SetData(data)
        pas.Abasis()
        #311.00403498262483
    """
    def __init__(self):
        self.data=[]
    def SetParam(self,param='n',val=0):
        """
        目的:paramで指定するパラメータの設定
        引数:
            param
                'n','mu','sigm','k'のいずれか
            val   設定する値
        """
        if param=='n':
            self.n=val
        elif param=='mu':
            self.mu=val
        elif param=='sigm':
            self.sigm=val
        elif param=='k':
            self.k=val
    def GetParam(self,param='n'):
        """
        目的:paramで指定するパラメータの呼びだし
        引数:
            param
                'n','mu','sigm','k'のいずれか
        """
        if param=='n':
            return self.n
        elif param=='mu':
            return self.mu
        elif param=='sigm':
            return self.sigm
        elif param=='k':
            return self.k
    def LowerLimit(self,mu=400,sigm=30,P=0.01):
        """
        目的:規定する破損確率Pに対する下限界許容応力
        """ 
        up=norm.ppf(1-P)
        s_a=mu-up*sigm
        return s_a
    def SetData(self,data):
        """
        目的:強度サンプルデータの読み込み
        """
        self.data=data
        self.n=len(data)
        self.mu=statistics.mean(data) #標本平均値
        self.sigm=statistics.stdev(data) #標本標準偏差
    def Basis1(self,P,gamma):
        """
        目的:許容値の計算(平均値，標準偏差ともに未知)
        """        
        up=norm.ppf(1-P)
        self.k=nct.ppf(1-gamma,self.n-1,np.sqrt(self.n)*up)/np.sqrt(self.n)
        return self.mu-self.k*self.sigm
    def Basis2(self,P,gamma,sigm):
        """
        目的:許容値の計算(平均値のみ未知)
        """        
        self.k=norm.ppf(1-P)+norm.ppf(1-gamma)/np.sqrt(self.n)
        s_a=self.mu-self.k*sigm
        return s_a
    def Abasis(self):
        """
        目的:A許容値の計算(平均値，標準偏差ともに未知)
        """
        P=0.01
        gamma=0.05
        s_a=self.Basis1(P,gamma)
        return s_a
    def Bbasis(self):
        """
        目的:B許容値の計算(平均値，標準偏差ともに未知)
        """
        P=0.1
        gamma=0.05
        s_a=self.Basis1(P,gamma)
        return s_a
    def Abasis2(self):
        """
        目的:A許容値の計算(平均値のみ未知)
        引数:
            sigm   既知の標準偏差
        """
        P=0.01
        gamma=0.05
        sigm=self.sigm
        s_a=self.Basis2(P,gamma,sigm)
        return s_a
    def Bbasis2(self):
        """
        目的:B許容値の計算(平均値のみ未知)
        引数:
            sigm   既知の標準偏差
        """
        P=0.1
        gamma=0.05
        sigm=self.sigm
        s_a=self.Basis2(P,gamma,sigm)
        return s_a
    def AS(self,spec='A1',sigm=20,P=0.1,gamma=0.05):
        """
        目的:
            汎用許容値評価
        引数:
            spec
                'A1':A許容値(平均値，標準偏差ともに未知)
                'A2':A許容値(平均値のみ未知,sigm指定)
                'B1':B許容値(平均値，標準偏差ともに未知)
                'B2':B許容値(平均値のみ未知,sigm指定)
                'C1':任意の許容値(平均値，標準偏差ともに未知,P,gamma指定)
                'C2':任意の許容値(平均値のみ未知,sigm,P,gamma指定)
            sigm 既知の標準偏差(A2,B2,C2の場合指定)
            P    破損確率
            gamma 1-信頼水準
        """
        if spec=='A1':
            P=0.01
            gamma=0.05
            s_a=self.Basis1(P,gamma)
        elif spec=='A2':
            P=0.01
            gamma=0.05
            s_a=self.Basis2(P,gamma,sigm)
        elif spec=='B1':
            P=0.1
            gamma=0.05
            s_a=self.Basis1(P,gamma)
        elif spec=='B2':
            P=0.1
            gamma=0.05
            s_a=self.Basis2(P,gamma,sigm)
        elif spec=='C1':
            s_a=self.Basis1(P,gamma)
        elif spec=='C2':
            s_a=self.Basis2(P,gamma,sigm)
        return s_a
class PAStressL:
    """    
    目的:
        許容応力の確率論的決定支援(対数正規分布)
    """
    def __init__(self):
        self.data=[]
        self.pas=PAStress()
    def SetData(self,data):
        """
        目的:強度サンプルデータの読み込み
        """
        self.data=np.log(data)
        self.pas.SetData(self.data)
    def Abasis(self):
        """
        目的:A許容値の計算(μL，σLともに未知)
        """
        x=self.pas.Abasis()
        return np.exp(x)
    def Bbasis(self):
        """
        目的:B許容値の計算(μL，σLともに未知)
        """
        x=self.pas.Bbasis()
        return np.exp(x)
    def Abasis2(self,sigm=2):
        """
        目的:A許容値の計算(μLのみ未知)
        引数:
            sigm: 母集団の対数データ標準偏差
        """
        x=self.pas.Abasis2(sigm)
        return np.exp(x)
    def Bbasis2(self,sigm=2):
        """
        目的:B許容値の計算(μLのみ未知)
        引数:
            sigm: 母集団の対数データ標準偏差
        """
        x=self.pas.Bbasis2(sigm)
        return np.exp(x)
from scipy.stats import chi2
class PAStressW:
    """    
    目的:
        許容応力の確率論的決定支援(形状母数既知の二母数ワイブル分布)
    """
    def __init__(self):
        self.data=[]
    def SetData(self,data,alpha):
        """
        目的:強度サンプルデータの読み込み
        data  データ配列
        alpha 既知の形状母数
        """
        self.alpha=alpha
        dd=np.array(data)
        self.data=dd
        self.beta=np.mean(dd**alpha)**(1./alpha)
        return self.beta
    def kAr(self,P,gam,n):
        """
        目的:k**alphaの計算
        　　　kは(kAr)**(1/self.alpha)により求まる
        """
        a=chi2.ppf(q=1-gam,df=2*n)
        return -a/(2*n*np.log(1-P))
    def Abasis(self):
        """
        目的:A許容値の計算
        """
        P=0.01
        gamma=0.05
        n=len(self.data)
        k=self.kAr(P,gamma,n)**(1./self.alpha)
        return self.beta/k
    def Bbasis(self):
        """
        目的:B許容値の計算
        """
        P=0.1
        gamma=0.05
        n=len(self.data)
        k=self.kAr(P,gamma,n)**(1./self.alpha)
        return self.beta/k