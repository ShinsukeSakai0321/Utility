"""
衝撃問題関連ユーティリティ　by S.Sakai
"""
from Utility import UnitConv as uc
import numpy as np
import math
import matplotlib.pyplot as plt
class WES:
    """
    溶接協会Wes式
    例題:三つのひずみ速度下での引張強さの計算(室温:降伏強さ228,引張強さ349)
        m1=SySt(228,349)
        St=[m1.St(1e-3),m1.St(1e-2),m1.St(1e-1)]
    """
    def __init__(self,Sy0,St0,T0=273+20,e0=0.2/(60),E=2.05E5):
        self.Sy=0.0
        self.Su=0.0
        self.Sy0=Sy0
        self.St0=St0
        self.T0=T0
        self.e0=e0
        self.E=E
    def Sy(self,e,T=273+20):
        c1=(1/(T*math.log(1e8/e))-1/(self.T0*math.log(1e8/self.e0)))
        val=self.Sy0*math.exp(8e-4*self.T0*(self.Sy0/self.E)**(-1.5)*c1)
        return val
    def St(self,e,T=273+20):
        c1=(1/(T*math.log(1e9/e))-1/(self.T0*math.log(1e9/self.e0)))
        val=self.St0*math.exp(8e-4*self.T0*(self.St0/self.E)**(-1.5)*c1)
        return val
class BRL:
    def eval(self,M,V,D,Ks=1.0):
        """
        目的:BRL式による貫通限界板厚の評価
        入力:
            M 飛来物質量[kg]
            V 衝突速度[m/s]
            Ks 鋼材のグレードによる係数
            D 飛来物有効直径[m]
        出力:
            T 貫通限界板厚
        
        """
        assert D!=0
        T=(0.5*M*V*V/(1.44e9*Ks*Ks*D*D*D)*D**1.5)**(1/1.5)
        return T
    def exam1(self):
        M=8.3
        K=0.84
        D=uc.Length().Conv('cm','m',10)
        V=80
        return M,K,D,V
class JohnsonCook:
    """
    Johnson-Cook則の計算
    登録材料 'DH36Steel''OFHC COPPER''ARMCO IRON''AISI-1045Steel''CARTRIDGE BRASS''NICKEL 200'
        'CARPENTER ELECTRICAL IRON''2024-T351 ALUMINUM''TUNGSTEN ALLOY''1006 STEEL''7039 ALUMINUM'
        '4340 STEEL''S-7 STEEL''DU-.75Ti''SM400C_1''SM400C_2'
    利用法:
        jc=JohnsonCook()
        jc.Material('DH36Steel')
        e_d0=0.01
        ep_d=0.2
        T=273+15
        s=jc.JC(ep,ep_d,T,e_d0=e_d0)
        print(s)
    """
    def __init__(self):
        self.A=0
        self.B=0
        self.n=0
        self.C=0
        self.m=0
        self.Tm=0
        self.T0=0
        self.ref=0
    def Reference(self):
        if self.ref==0:
            ref="G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'"
        elif self.ref==1:
            ref="M.Murugesan,D.W.Jung,'Johnson Cook Material and Failure Model Parameters Estimation of AISI-1045 Medium Carbon Steel for Metal Forming Applications',Materials 2019, 12, 609; doi:10.3390/ma12040609"
        elif self.ref==2:
            ref="HPI SIL3"
        elif self.ref==3:
            ref="HPI SIL2"
        return ref
    def Material(self,mat):
        assert mat!=''
        if mat=='SM400C_1':
            self.A=3.101e2
            self.B=5.694e2
            self.n=4.931e-1
            self.C=4.8573e-4
            self.e_d0=1e-3
            self.ref=3
        if mat=='SM400C_2':
            self.A=3.101e2
            self.B=5.289e2
            self.n=4.767e-1
            self.C=3.575e-4
            self.e_d0=1e-1
            self.ref=3
        if mat=='DH36Steel':
            self.A=1020
            self.B=1530
            self.n=0.4
            #self.e0_d=0.1
            self.C=0.015
            self.m=0.32
            #self.T0=50
            self.Tm=1773
            self.ref=2
        if mat=='OFHC COPPER':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=90
            self.B=292
            self.n=0.31
            self.C=0.025
            self.m=1.09
            self.Tm=1356
            self.T0=273+15
            self.ref=0
        if mat=='ARMCO IRON':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=175
            self.B=380
            self.n=0.32
            self.C=0.06
            self.m=0.55
            self.Tm=1811
            self.T0=273+15
            self.ref=0
        if mat=='CARTRIDGE BRASS':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=112
            self.B=505
            self.n=0.42
            self.C=0.009
            self.m=1.68
            self.Tm=1189
            self.T0=273+15
            self.ref=0
        if mat=='NICKEL 200':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=163
            self.B=648
            self.n=0.33
            self.C=0.006
            self.m=1.44
            self.Tm=1726
            self.T0=273+15
            self.ref=0
        if mat=='CARPENTER ELECTRICAL IRON':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=290
            self.B=339
            self.n=0.4
            self.C=0.055
            self.m=0.55
            self.Tm=1811
            self.T0=273+15
            self.ref=0
        if mat=='1006 STEEL':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=350
            self.B=275
            self.n=0.36
            self.C=0.022
            self.m=1.0
            self.Tm=1811
            self.T0=273+15
            self.ref=0
        if mat=='2024-T351 ALUMINUM':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=265
            self.B=426
            self.n=0.34
            self.C=0.015
            self.m=1.00
            self.Tm=775
            self.T0=273+15
            self.ref=0
        if mat=='7039 ALUMINUM':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=337
            self.B=343
            self.n=0.41
            self.C=0.01
            self.m=1.00
            self.Tm=877
            self.T0=273+15
            self.ref=0
        if mat=='4340 STEEL':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=792
            self.B=510
            self.n=0.26
            self.C=0.014
            self.m=1.03
            self.Tm=1793
            self.T0=273+15
            self.ref=0
        if mat=='S-7 STEEL':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=1539
            self.B=477
            self.n=0.18
            self.C=0.012
            self.m=1.0
            self.Tm=1763
            self.T0=273+15
            self.ref=0
        if mat=='TUNGSTEN ALLOY':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=1506
            self.B=177
            self.n=0.12
            self.C=0.016
            self.m=1.00
            self.Tm=1723
            self.T0=273+15
            self.ref=0
        if mat=='DU-.75Ti':
            #G.R.Johnson,W.H.Cook,'A Constitutive Medel And Data for Metals Subjected to Large Strains, High Strain Rates and Heigh Temeratures'
            self.A=1079
            self.B=1120
            self.n=0.25
            self.C=0.007
            self.m=1.00
            self.Tm=1473
            self.T0=273+15
            self.ref=0
        if mat=='AISI-1045Steel':
            #M.Murugesan,D.W.Jung,'Johnson Cook Material and Failure Model Parameters Estimation of AISI-1045 Medium Carbon Steel for Metal Forming Applications'Materials 2019, 12, 609; doi:10.3390/ma12040609
            self.A=50.103
            self.B=176.091
            self.n=0.5176
            self.C=0.095
            self.m=0.6622
            self.Tm=1623
            self.T0=1223
            self.ref=1
    def JC(self,e,e_d,T=0,e_d0=1):
        if T==0:
            th=0.0
        else:
            th=(T-self.T0)/(self.Tm-self.T0)
        c1=(self.A+self.B*e**self.n)
        c2=(1+self.C*math.log(e_d/e_d0))
        if T==0:
            c3=1.0
        else:
            c3=(1-th**self.m)
        return c1*c2*c3