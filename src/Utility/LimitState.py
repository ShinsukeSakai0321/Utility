"""
    限界状態関数法サポートモジュール
　　Copyright 2021, S.Sakai All rights reserved
"""
import numpy as np
from scipy.stats import norm
class Dbase:
    """
    目的:任意分布管理のための基底クラス
    """
    def __init__(self,mu,sigmma):
        self.muX=mu
        self.sigmmaX=sigmma
    def GetMeanSig(self):
        return [self.muX,self.sigmmaX]
    def DushCalc(self,X):
        self.sigmmaDush=norm.pdf(norm.ppf(self.CDF(X))/self.PDF(X))
        self.muDush=X-(norm.ppf(self.CDF(X)))*self.sigmmaDush
    def CDF(self,X):
        return norm.cdf(X)
    def PDF(self,X):
        return norm.pdf(X)
    def Eq(self,X):
        print('Specify Eq()!')
import math
from scipy.stats import lognorm
from scipy.stats import gumbel_r
from scipy.stats import weibull_min
class GNormal(Dbase):
    def __init__(self,mu,sigmmaX):
        super().__init__(mu,sigmmaX)
    def Eq(self,X):
        return [self.muX,self.sigmmaX]
class GLognormal(Dbase):
    """
       zeta:σ
       lambd:scale
    """
    def __init__(self,mu,sigmmaX):
        super().__init__(mu,sigmmaX)
    def Eq(self,X):
        zeta=math.sqrt(math.log(1+(self.sigmmaX/self.muX)**2))
        lambd=math.log(self.muX)-0.5*zeta*zeta
        scale=math.exp(lambd)
        phi=norm.pdf(norm.ppf(lognorm(s=zeta,scale=scale).cdf(X)))
        fxi=lognorm(s=zeta,scale=scale).pdf(X)
        sigm=phi/fxi
        mu=X-norm.ppf(lognorm(s=zeta,scale=scale).cdf(X))*sigm
        return [mu,sigm]
    def Param(self):
        self.zeta=math.sqrt(math.log(1+(self.sigmmaX/self.muX)**2))
        self.lambd=math.log(self.muX)-0.5*self.zeta*self.zeta
        return [self.zeta,self.lambd]
    def SetParam(self,lambd,zeta):
        self.muX=math.exp(lambd+zeta*zeta/2.0)
        self.sigmmaX=math.sqrt(self.muX*self.muX*(math.exp(eta*zeta)-1.0))
        return [self.muX,self.sigmmaX]
class GGumbel(Dbase):
    def __init__(self,mu,sigmmaX):
        super().__init__(mu,sigmmaX)
    def Eq(self,X):
        eta=math.sqrt(6.)*self.sigmmaX/math.pi #shape
        mu=self.muX-0.57722*eta #loc
        phi=norm.pdf(norm.ppf(gumbel_r.cdf(X,mu,eta)))
        fXi=gumbel_r.pdf(X,mu,eta)
        sigm=phi/fXi
        mu=X-norm.ppf(gumbel_r.cdf(X,mu,eta)) * sigm
        return [mu,sigm]
    def Param(self):
        eta=math.sqrt(6.)*self.sigmmaX/math.pi #shape
        mu=self.muX-0.57722*eta #loc
        return [mu,eta]
    def SetParam(self,mu,eta):
        self.sigmmaX=math.pi/math.sqrt(6.)*eta
        self.muX=mu+0.57722*eta
        return [self.muX,self.sigmmaX]
    def MuSigm(self,mu,eta):
        return [mu+0.57722*eta, math.pi*eta/math.sqrt(6.)]
class GWeibull(Dbase):
    def __init__(self,mu,sigmmaX):
        super().__init__(mu,sigmmaX)
    def Eq(self,X):
        self.FactorCalc()
        super().DushCalc(X)
        return [self.muDush,self.sigmmaDush]
    def FactorCalc(self):
        cov=self.muX/self.sigmmaX
        self.ar=cov**(-1.08) #Prof.Ichikawa's theory
        self.beta=self.muX/math.gamma(1.0+1/self.ar)
    def PDF(self,X):
        return weibull_min.pdf(X, self.ar, scale=self.beta)
    def CDF(selff,X):
        return weibull_min.cdf(X, self.ar, scale=self.beta)
class Lbase:
    """
    目的:限界状態関数定義のための基底クラス
    """
    def __init__(self,n):
        self.n=n
        self.X=np.zeros(n)
        self.dGdX=np.zeros(n)
        self.g=0
    def GetN(self):
        return self.n
    def GetX(self):
        return self.X
    def SetX(self,X):
        self.X=X
    def GetG(self):
        return self.g
    def SetG(self,g):
        self.g=g
    def GetdGdX(self):
        return self.dGdX
    def SetdGdX(self,x):
        self.dGdX=x
from sympy import sympify,Symbol,factor,diff,pprint,expand
class GeneralG(Lbase):
    def __init__(self,gg,var):
        self.n=len(var)
        super().__init__(self.n)
        self.gg=gg
        self.var=var
    def gcalc(self):
        expr=self.gg
        X=super().GetX()
        for i in range(self.n):
            str1=self.var[i]+'=X['
            str1=str1+str(i)+']'
            #eval(str1)
            exec(str1)
        super().SetG(eval(expr))
    def dGdXcalc(self):
        X=super().GetX()
        dGdX=super().GetdGdX()
        for i in range(self.n):
            str1=self.var[i]+'=X['
            str1=str1+str(i)+']'
            exec(str1)
        expr=sympify(self.gg,locals={'S': Symbol('S')})
        for i in range(self.n):
            dstr1=diff(expr,self.var[i])
            dGdX[i]=eval(str(dstr1))
        super().SetdGdX(dGdX) 
import pandas as pd
class LSFM:
    """
    目的:信頼性解析管理のための基底クラス
    method一覧
    RF()      Rackvitz Fiessler法による設計点探索
    GetBeta() 信頼性指標βの取得
    GetAlpha()感度ベクトルの取得
    GetPOF()  破損確率の取得
    GetDP()   設計点の取得
    GetConv() RF法の収束回数の取得
    GetPSF()  部分安全係数の取得
    """
    def __init__(self,n,Mu,sigmmaX,dist):
        self.n=n
        self.muX=Mu
        self.sigmmaX=sigmmaX
        self.dist=dist
        self.alphai=np.zeros(n)
        self.beta=0.0
        self.Distr=[]
        #self.lim=Lbase
        for i in range(n):
            if self.dist[i]=='normal':
                self.Distr.append(GNormal(self.muX[i],self.sigmmaX[i]))
            if self.dist[i]=='lognormal':
                self.Distr.append(GLognormal(self.muX[i],self.sigmmaX[i]))
            if self.dist[i]=='gumbel':
                self.Distr.append(GGumbel(self.muX[i],self.sigmmaX[i]))
    def GetN(self):
        return self.n
    def GetMu(self):
        return self.muX
    def SetMu(self,aa):
        self.muX=aa
    def RF(self,Xstart):
        """
        Racwitz-Fiessler algoritm
        Xstart: list of starting point
        """
        betaold=40
        delta=1e-6
        munormX=np.zeros(self.n)
        sigmmanormX=np.zeros(self.n)
        #X=self.muX.copy()
        X=Xstart
        for i in range(100):
            for j in range(self.n):
                Valu=self.Distr[j].Eq(X[j])
                munormX[j]=Valu[0]#非正規分布の正規化後の平均値
                sigmmanormX[j]=Valu[1]#非正規分布の正規化後の標準偏差
            Xdush=(X-munormX)/sigmmanormX
            Xdush_old=Xdush
            self.lim.SetX(X)
            self.lim.gcalc()
            self.lim.dGdXcalc()
            g=self.lim.GetG()
            dgdX=self.lim.GetdGdX()
            dgdXdush=dgdX*sigmmanormX
            A = 1 / sum(dgdXdush * dgdXdush) * (sum(dgdXdush * Xdush) - g)
            Xdush = A * dgdXdush
            self.alphai = dgdXdush / math.sqrt(sum(dgdXdush * dgdXdush))
            Xdush_new = Xdush
            betanew = math.sqrt(sum(Xdush * Xdush))
            hantei = math.isnan(betanew)
            if hantei:
                betaold = betaold
            else:
                X = munormX + sigmmanormX * Xdush
                if abs(betaold - betanew) < delta:
                    break
                betaold=betanew
            X_t = munormX + sigmmanormX * Xdush_old
            self.lim.SetX(X_t)
            self.lim.gcalc()
            g_hantei = self.lim.GetG()
            deltan = 50
            dXdush = (Xdush_new - Xdush_old) / deltan
            for i1 in range(deltan):
                Xdush_t = Xdush_old + i1 * dXdush
                X_t = munormX + sigmmanormX * Xdush_t
                self.lim.SetX(X_t)
                self.lim.gcalc()
                g_hantein = self.lim.GetG()
                if math.isnan(g_hantei): #<----------------
                    return
                if g_hantei**2 > g_hantein**2:
                    Xdush = Xdush_t
                    g_hantei = g_hantein
                else:
                    g_hantei =g_hantei
            X = munormX + sigmmanormX * Xdush
        self.beta=betanew
        self.POF=norm.sf(betanew)
        self.DP=X
        self.ncon=i
    def GetSigm(self):
        return self.sigmmaX
    def GetBeta(self):
        return self.beta
    def GetAlpha(self):
        return self.alphai
    def GetPOF(self):
        return self.POF
    def GetDP(self):
        return self.DP
    def GetConv(self):
        return self.ncon
    def GetPSF(self):
        return self.DP/self.muX
    def DefineG(self,glim):
        self.lim=glim
    def GetLim(self):
        return self.lim
    def GetG(self):
        X=self.muX
        self.lim.SetX(X)
        self.lim.gcalc()
        return self.lim.GetG()
    def Gcalc(self,X):
        self.lim.SetX(X)
        self.lim.gcalc()
        return self.lim.GetG()
    def GetdGdX(self):
        return self.lim.GetdGdX()
class GeneralTreat(LSFM):
    """
    目的:限界状態関数定義を文字列で与える場合のGを定義するクラス
    """
    def __init__(self,g,var,dist,muX,sigmmaX):
        n=len(var)
        self.g=g
        self.var=var
        super().__init__(n,muX,sigmmaX,dist)
    def calc(self):
        gg=GeneralG(self.g,self.var)
        super().DefineG(gg)
        super().RF()
def reliability(g='r-s',var=['r','s'],dist=['normal','normal'],muX=[200,100],sigmmaX=[10,20]):
    """
    目的:限界状態関数定義を文字列で与える場合の信頼性解析を管理するクラス
    """
    aa=GeneralTreat(g,var,dist,muX,sigmmaX)
    aa.calc()
    return aa
class Gmanage(LSFM):
    """ユーザ定義のGをLSFMによる解析に接続する
    　　使い方
      　brl=Gmanage(n,muX,sigmmaX,dist,G)
        brl.RF()
        print('beta=',brl.GetBeta())
        print('alpha=',brl.GetAlpha())
        Gの部分にユーザ定義のクラス名を記述
    """
    def __init__(self,n,Mu,sigmmaX,dist,g):
        self.n=n
        gg=g(self.n)
        super().DefineG(gg)
        super().__init__(n,Mu,sigmmaX,dist)
class RelBase:
    """汎用信頼性評価のための基底クラス
    """
    def SetData(self,data):
        key=self.variable
        self.muX=[]
        self.cov=[]
        self.dist=[]
        """
        for aa in key:
            self.muX.append(data[aa]['mean'])
            self.cov.append(data[aa]['cov'])
            self.dist.append(data[aa]['dist'])
        self.sigmmaX = list(np.array(self.cov)*np.array(self.muX))
        """
        self.sigmmaX=[]
        for aa in key:
            self.dist.append(data[aa]['dist'])
            muX=data[aa]['mean']
            self.muX.append(muX)
            if 'cov' in data[aa]:
                sX=muX*data[aa]['cov']
            else:
                sX=data[aa]['sd']
            self.sigmmaX.append(sX)       
    def Reliability(self,data,start):
        n=len(self.variable)
        self.SetData(data)
        self.lsfm=Gmanage(n,self.muX,self.sigmmaX,self.dist,self.G)
        self.lsfm.RF(start)
    def GetDP(self):
        return self.lsfm.GetDP()
    def Gcalc(self,X):
        return self.lsfm.Gcalc(X)
    def Gmean(self):
        return self.lsfm.Gcalc(self.muX)
    def GetBeta(self):
        return self.lsfm.GetBeta()
    def GetAlpha(self):
        return self.lsfm.GetAlpha()
    def GetPOF(self):
        return self.lsfm.GetPOF()
    def GetVar(self):
        return self.variable
    def GetTitle(self):
        return self.title
from math import exp,log,sqrt
class GMetalLoss(Lbase):
    """
    目的:API法による局部減肉信頼性評価のための限界状態関数定義
    """
    def __init__(self,n):
        self.n=n
        super().__init__(self.n)
    def Settime(self,t):
        self.ttime=t
    def gcalc(self):
        X=super().GetX()
        ttime=self.ttime
        Cv=X[0]
        Papplied=X[1]
        tr=X[2]
        no=X[3]
        D=X[4]
        sigmmau=X[5]
        s=X[6]
        tmm=X[7]
        tc=tr-Cv*ttime
        lambd=1.2850*s/sqrt(D*tc)
        Mt=1.0010 - 0.014195 * lambd + 0.29090 * lambd**2 - 0.096420 * lambd**3 + 0.020890 * lambd**4 - 0.0030540 * lambd**5 + 2.9570e-4 * lambd**6 - 1.8462e-5 * lambd**7 + 7.1553e-7 * lambd**8 - 1.531e-8 * lambd**9 + 1.4656e-10 * lambd**10
        Rt =max(0, (tmm - Cv * ttime) / tc)
        RSF=Rt/(1.0-(1.0-Rt)/Mt)
        Pbi = (exp(no)/no**no) * (0.25 / (no + 0.227)) *log(1.0 + 2.0 * tc / D) * sigmmau
        g=RSF*Pbi-Papplied
        super().SetG(g)
    def dGdXcalc(self):
        X=super().GetX()
        dGdX=super().GetdGdX()
        ttime=self.ttime
        Cv=X[0]
        Papplied=X[1]
        tr=X[2]
        no=X[3]
        D=X[4]
        sigmmau=X[5]
        s=X[6]
        tmm=X[7]
        dGdX[0] = -0.5*exp(no)*sigmmau*ttime*(tmm-Cv*ttime)/(no**no*(no+0.227)*(tr-Cv*ttime)*(2.0*(tr-Cv*ttime)/D+1.0)*D*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*ttime*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))+0.25*exp(no)*sigmmau*ttime*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)**2*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)*((1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))*(-0.0091202875*s*ttime*D*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.30687942319875*s**3*ttime*D*((tr-Cv*ttime)*D)**((-5.0)/2.0)-0.0267500550265318*s**5*ttime*D*((tr-Cv*ttime)*D)**((-7.0)/2.0)-3.738252949569017E-4*s**7*ttime*D*((tr-Cv*ttime)*D)**((-9.0)/2.0)-6.581363123987437E-7*s**9*ttime*D*((tr-Cv*ttime)*D)**((-11.0)/2.0)+0.4803413525*s**2*ttime/((tr-Cv*ttime)**2*D)+0.113915008346112*s**4*ttime/((tr-Cv*ttime)**3*D**2)+0.00399384628042394*s**6*ttime/((tr-Cv*ttime)**4*D**3)+2.1277120825241515E-5*s**8*ttime/((tr-Cv*ttime)**5*D**4)+8.995322480552324E-9*s**10*ttime/((tr-Cv*ttime)**6*D**5))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)**2-(ttime/(tr-Cv*ttime)-ttime*(tmm-Cv*ttime)/(tr-Cv*ttime)**2)/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001))/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001))**2)
        dGdX[1] = -1.0    
        dGdX[2] = 0.5*exp(no)*sigmmau*(tmm-Cv*ttime)/(no**no*(no+0.227)*(tr-Cv*ttime)*(2.0*(tr-Cv*ttime)/D+1.0)*D*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)**2*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)*((1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))*(0.0091202875*s*D*((tr-Cv*ttime)*D)**((-3.0)/2.0)+0.30687942319875*s**3*D*((tr-Cv*ttime)*D)**((-5.0)/2.0)+0.0267500550265318*s**5*D*((tr-Cv*ttime)*D)**((-7.0)/2.0)+3.738252949569017E-4*s**7*D*((tr-Cv*ttime)*D)**((-9.0)/2.0)+6.581363123987437E-7*s**9*D*((tr-Cv*ttime)*D)**((-11.0)/2.0)-0.4803413525*s**2/((tr-Cv*ttime)**2*D)-0.113915008346112*s**4/((tr-Cv*ttime)**3*D**2)-0.00399384628042394*s**6/((tr-Cv*ttime)**4*D**3)-2.1277120825241515E-5*s**8/((tr-Cv*ttime)**5*D**4)-8.995322480552324E-9*s**10/((tr-Cv*ttime)**6*D**5))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)**2-(tmm-Cv*ttime)/((tr-Cv*ttime)**2*(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001))**2)
        dGdX[3] = 0.25*exp(no)*(-log(no)-1)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))+0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)**2*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))
        dGdX[4] = -0.5*exp(no)*sigmmau*(tmm-Cv*ttime)/(no**no*(no+0.227)*(2.0*(tr-Cv*ttime)/D+1.0)*D**2*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))*log(2.0*(tr-Cv*ttime)/D+1.0)*(0.0091202875*s*(tr-Cv*ttime)*((tr-Cv*ttime)*D)**((-3.0)/2.0)+0.30687942319875*s**3*(tr-Cv*ttime)*((tr-Cv*ttime)*D)**((-5.0)/2.0)+0.0267500550265318*s**5*(tr-Cv*ttime)*((tr-Cv*ttime)*D)**((-7.0)/2.0)+3.738252949569017E-4*s**7*(tr-Cv*ttime)*((tr-Cv*ttime)*D)**((-9.0)/2.0)+6.581363123987437E-7*s**9*(tr-Cv*ttime)*((tr-Cv*ttime)*D)**((-11.0)/2.0)-0.4803413525*s**2/((tr-Cv*ttime)*D**2)-0.113915008346112*s**4/((tr-Cv*ttime)**2*D**3)-0.00399384628042394*s**6/((tr-Cv*ttime)**3*D**4)-2.1277120825241515E-5*s**8/((tr-Cv*ttime)**4*D**5)-8.995322480552324E-9*s**10/((tr-Cv*ttime)**5*D**6))/(no**no*(no+0.227)*(tr-Cv*ttime)*(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)**2*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001))**2)
        dGdX[5] = 0.25*exp(no)*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.0107000220106127*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.0569575041730562*s**4/((tr-Cv*ttime)**2*D**2)+0.00133128209347465*s**6/((tr-Cv*ttime)**3*D**3)+5.319280206310379E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))
        dGdX[6] = -0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))*log(2.0*(tr-Cv*ttime)/D+1.0)*(-0.018240575/sqrt((tr-Cv*ttime)*D)-0.6137588463975*s**2*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.053500110053064*s**4*((tr-Cv*ttime)*D)**((-5.0)/2.0)-7.4765058991380345E-4*s**6*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.3162726247974875E-6*s**8*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.960682705*s/((tr-Cv*ttime)*D)+0.22783001669222*s**3/((tr-Cv*ttime)**2*D**2)+0.0079876925608479*s**5/((tr-Cv*ttime)**3*D**3)+4.255424165048303E-5*s**7/((tr-Cv*ttime)**4*D**4)+1.7990644961104649E-8*s**9/((tr-Cv*ttime)**5*D**5))/(no**no*(no+0.227)*(tr-Cv*ttime)*(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.010700022010613*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.056957504173056*s**4/((tr-Cv*ttime)**2*D**2)+0.0013312820934746*s**6/((tr-Cv*ttime)**3*D**3)+5.3192802063103787E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)**2*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.010700022010613*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.056957504173056*s**4/((tr-Cv*ttime)**2*D**2)+0.0013312820934746*s**6/((tr-Cv*ttime)**3*D**3)+5.3192802063103787E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001))**2)
        dGdX[7] = 0.25*exp(no)*sigmmau*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.010700022010613*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.056957504173056*s**4/((tr-Cv*ttime)**2*D**2)+0.0013312820934746*s**6/((tr-Cv*ttime)**3*D**3)+5.3192802063103787E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)))-0.25*exp(no)*sigmmau*(tmm-Cv*ttime)*log(2.0*(tr-Cv*ttime)/D+1.0)/(no**no*(no+0.227)*(tr-Cv*ttime)**2*(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.010700022010613*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.056957504173056*s**4/((tr-Cv*ttime)**2*D**2)+0.0013312820934746*s**6/((tr-Cv*ttime)**3*D**3)+5.3192802063103787E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001)*(1.0-(1.0-(tmm-Cv*ttime)/(tr-Cv*ttime))/(-0.018240575*s/sqrt((tr-Cv*ttime)*D)-0.2045862821325*s**3*((tr-Cv*ttime)*D)**((-3.0)/2.0)-0.010700022010613*s**5*((tr-Cv*ttime)*D)**((-5.0)/2.0)-1.0680722713054335E-4*s**7*((tr-Cv*ttime)*D)**((-7.0)/2.0)-1.462525138663875E-7*s**9*((tr-Cv*ttime)*D)**((-9.0)/2.0)+0.4803413525*s**2/((tr-Cv*ttime)*D)+0.056957504173056*s**4/((tr-Cv*ttime)**2*D**2)+0.0013312820934746*s**6/((tr-Cv*ttime)**3*D**3)+5.3192802063103787E-6*s**8/((tr-Cv*ttime)**4*D**4)+1.7990644961104648E-9*s**10/((tr-Cv*ttime)**5*D**5)+1.001))**2)
        super().SetdGdX(dGdX)
class MetalLoss(LSFM):
    """
    目的:API法による局部減肉信頼性評価のための管理クラス
    """
    def __init__(self,n,Mu,sigmmaX,dist,ttime):
        self.n=n
        self.ttime=ttime
        gg=GMetalLoss(self.n)
        gg.Settime(self.ttime)
        super().DefineG(gg)
        super().__init__(n,Mu,sigmmaX,dist)
    def calc(self):
        super().RF()
    def PSFCalc(self):
        psf=self.muX/super().GetDP()
        psf[0]=1/psf[0]
        psf[1]=1/psf[1]
        return psf
    def G(self):
        return super().GetG()  
from pyDOE import *
class LHSbase:
    """
    目的:LHS+USの乱数発生を管理するための基底クラス
    RSモデルの例
        from Utility import LimitState as ls
        class rnd_RS(ls.LHSbase):
            def __init__(self,nv):
                super().__init__(nv)
            def g(self,rnd):  #限界状態関数gを定義する。rnd:発生された乱数列。変数の順番に格納されている。
                    rr=rnd[:,0]  #変数の順番に乱数列を取り出していく
                    ss=rnd[:,1]
                    gval=rr-ss  #限界状態関数の計算結果リスト
                    return gval
        #### プログラム例 ####
        data={"r":{"mean":170,"cov":20/170,"dist":"normal"},
             "s":{"mean":100,"cov":20/100,"dist":"normal"},
             }
        nv=2 #変数の数
        n=300 #サンプル点数
        k=3 #サンプリングの領域の広さ
        lhb=rnd_RS(nv) #インスタンスの生成
        lhb.SetData(data) #データセット
        rnd,t=lhb.Calc(n,k) #乱数発生
        X_std=(rnd-lhb.Means())/lhb.Stdvs() #データの基準化
        print(lhb.gMean()) #平均値でのg値の値

    """
    def __init__(self,nv):
        self.nv=nv
    def SetData(self,data):
        self.data=data
        means=[]
        stdvs=[]
        for key in data:
            mean=data[key]['mean']
            std=data[key]['cov']*mean
            means.append(mean)
            stdvs.append(std)
        self.means=means
        self.stdvs=stdvs
    def Means(self):
        return self.means
    def Stdvs(self):
        return self.stdvs
    def gMean(self):
        num=self.nv
        val=np.zeros((num,2))
        for i in range(num):
            val[0,i]=self.means[i]
        return self.g(val)[0]
    def Calc(self,n,k):
        design = lhs(self.nv, samples=n,criterion='maximin')
        u=design
        for i in range(self.nv):
            u[:, i] = k*self.stdvs[i]*(2*design[:,i]-1)+self.means[i]
        rnd=u
        t=(np.sign(self.g(rnd))+1)/2
        return rnd,t             