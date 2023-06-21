from Utility import LimitState as ls
from sympy import *
import numpy as np
class penMed(ls.RelBase):
    """
    目的:貫通評価クラスが継承する直前のクラス
    """
    def MakeContour(self,data,cdata):
        """
        目的:JSONデータdata内の二つのパラメータに関するPOF等高線データ作成
            data:   貫通評価計算のための入力データ(JSON形式)
            cdata:  計算する格子点に関する情報を格納するJSONデータ
                例: 変数bとmに関する等高線データを発生するとき
                    cdata={'b':{'min':0.006,'max':0.012,'div':100},
                        'm':{'min':1,'max':3,'div':100}}
        戻り値: X,Y,Z
            plt.pcolormesh(X, Y, Z[i], cmap='hsv')などにより等高線描画する
            ここで
                i=0: PoF
                i=1: Beta
                i=2; Alpha_0
                i=3: Alpha_1
                .
                .
            Alpha_iは感度値を示し、i=0,1,...はsuper().GetVar()により出力される変数の順番と対応している
        """
        key=list(cdata.keys())
        x=np.arange(cdata[key[0]]['min'],cdata[key[0]]['max'],(cdata[key[0]]['max']-cdata[key[0]]['min'])/cdata[key[0]]['div'])
        y=np.arange(cdata[key[1]]['min'],cdata[key[1]]['max'],(cdata[key[1]]['max']-cdata[key[1]]['min'])/cdata[key[1]]['div'])
        X, Y = np.meshgrid(x, y)
        #ZZ=[]
        ZZ=[[] for i in range(len(super().GetVar())+2)]
        for iy in range(len(y)):    
            yy=Y[:,0][iy]
            data[key[1]]['mean']=yy
            zPoF=[]
            zBeta=[]
            zAlpha=[[] for i in range(len(super().GetVar()))]
            for ix in range(len(x)):
                xx=X[0][ix]
                data[key[0]]['mean']=xx
                super().Reliability(data)
                zPoF.append(super().GetPOF())
                zBeta.append(super().GetBeta())
                for ii in range(len(super().GetVar())):
                    zAlpha[ii].append(super().GetAlpha()[ii])
            ZZ[0].append(zPoF)
            ZZ[1].append(zBeta)
            for ii in range(len(super().GetVar())):
                ZZ[ii+2].append(zAlpha[ii])
        return X,Y,ZZ

        


    def SaveRange(self,aa):
        """
        目的:適用範囲データの保存
        """
        self.Range=aa
    def ShowRange(self):
        """
        目的:適用範囲データの表示
        """
        return self.Range
    def check(self,cond,val):
        """
        目的:適用範囲内であるかどうかのチェック
        　　使い方:  継承したクラスから
          　　　super().check('b/d',b/d)
             など
        """
        min_r=self.Range[cond][0]
        max_r=self.Range[cond][1]
        if val >= min_r and val<=max_r:
            print('**Validation of [',cond,'] satisfied**')
            return
        print('**Validation of [',cond,'] not satisfied**:',',Value=',val)
    def GetMean(self,data):
        """
        目的:Jsonデータ　dataから、変数の平均値を取り出しリストを戻す
        """
        dmean=[]
        for aa in data.keys():
            if 'cov' in data[aa].keys() or 'sd' in data[aa].keys():
                dmean.append(data[aa]['mean'])
        return dmean



"""
貫通評価モジュール
"""
################################
#             BRL              #
################################
class BRL(penMed):
    """
    Ballistic Research Laboratories (BRL) model (1968)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    ***constants***
    Limp    length of impactor
    Lsh     unsupported shield panel span
    Su      ultimate tensile strength of shield material
    """
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='BRL Formula'
        val_range={
            'v_bl':[57,270],
            'Limp/d':[1.25,8],
            'b/d':[0.1,1.0],
            'Lsh/d':[8,35],
            'Su':[315,500]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Limp=data['Limp']['mean']
        Lsh=data['Lsh']['mean']
        Su=data['Su']['mean']
        a7=5.37
        v_bl=a7*1e4*(b*d)**0.75/m**0.5
        super().check('v_bl',v_bl)
        super().check('Limp/d',Limp/d)
        super().check('b/d',b/d)
        super().check('Lsh/d',Lsh/d)
        super().check('Su',Su)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
            b,d,m,v=symbols('b d m v')
            a7=5.37
            var('a7')
            g=a7*1e4*(b*d)**0.75/m**0.5-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            a7=5.37
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            a7=5.37
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            super().SetdGdX(dGdX)
################################
#             DeMarre          #
################################
class DeMarre(penMed):
    """
    De Marre formula (Herrmann and Jones,1961)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='De Marre Formula'
        val_range={
            'v_bl':[200,900],
            'm':[0.1,50]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        v_bl=0.4311e5*d**0.75*b**0.7/m**0.5
        super().check('v_bl',v_bl)
        super().check('m',m)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
            b,d,m,v=symbols('b d m v')
            g=0.431e5*d**0.75*b**0.7/m**0.5-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            super().SetdGdX(dGdX)
################################
#             THOR             #
################################
import numpy as np

class THOR(penMed):
    """
    THOR equation (Crull and Swisdak, 2005)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    th  angle between a normal vector to a shield surface and the direction of impactor
    """
    C1=0
    a1=0
    b1=0
    g1=0
    def __init__(self):
        self.variable=['b','d','m','v','th']
        self.title='THOR Equations'
        self.tab={"Magnesium":{"C1":6.349,"a1":1.004,"b1":-1.076,"g1":0.966},
             "Aluminum":{"C1":6.185,"a1":0.903,"b1":-0.941,"g1":1.098},
             "CastIron":{"C1":10.153,"a1":2.186,"b1":-2.204,"g1":2.156},
             "Titanium":{"C1":7.552,"a1":1.325,"b1":-1.314,"g1":1.643},
             "FaceSteel":{"C1":7.694,"a1":1.191,"b1":-1.397,"g1":1.747}, 
             "MildSteel":{"C1":6.523,"a1":0.906,"b1":-0.963,"g1":1.286}, 
             "HardSteel":{"C1":6.601,"a1":0.906,"b1":-0.963,"g1":1.286},
             "Copper":{"C1":14.065,"a1":3.476,"b1":-3.687,"g1":4.27},
             "Lead":{"C1":10.955,"a1":2.735,"b1":-2.753,"g1":3.59}
            }
        super().SaveRange('Validation process is not defined.')
        super().SaveVariable(self.variable)
    def MatList(self):
        """
        目的:登録されている材料名リストを返す
        """
        return list(self.tab.keys())
    def setMaterial(self,mat):
        """
        目的:材料の設定
           mat: "Magnesium","Aluminum","CastIron","Titanium","FaceSteel","MildSteel","HardSteel","Copper","Lead"
        """
        global C1,a1,b1,g1
        C1=self.tab[mat]["C1"]
        a1=self.tab[mat]["a1"]
        b1=self.tab[mat]["b1"]
        g1=self.tab[mat]["g1"]
    class G(ls.Lbase):
        def __init__(self,n):
            global C1,a1,b1,g1
            self.n=n
            super().__init__(self.n)
            b,d,m,v,th=symbols('b d m v th')
            A=np.pi*d*d/4
            var('A C1 a1 b1 g1')
            g=0.3048*10**C1*(61024*b*A)**a1*(15432.4*m)**b1*(1/cos(th))**g1-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,v))
            self.d4=str(diff(g,th))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            th=X[4]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            th=X[4]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            super().SetdGdX(dGdX)
################################
#             Ohte             #
################################
class Ohte(penMed):
    """
    Ohte et al. Formula (Ohte et al., 1982)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='Ohte et al. Formula'
        val_range={
            'v_bl':[25,180],
            'm':[3,50],
            'Su':[490,637],
            'b':[7,38],
            'Lsh/b':[39,1e4],
            'd':[39,1e4]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        v_bl=7.67e4*(b*d)**0.75/m**0.5
        super().check('v_bl',v_bl)
        super().check('m',m)
        super().check('Su',Su)
        super().check('b',b)
        super().check('Lsh/b',Lsh/b)
        super().check('d',d)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
            b,d,m,v=symbols('b d m v')
            g=7.67e4*(b*d)**0.75/m**0.5-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            super().SetdGdX(dGdX)
################################
#             SRI              #
################################
class SRI(penMed):
    """
    Stanford Research Institute (SRI) correlation (1963)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    Lsh unsupported shield panel span
    Su  ultimate tensile strength of shield material
    Limp    length of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','v','Lsh','Su']
        self.title='Ohte et al. Formula'
        val_range={
            'v_bl':[21,122],
            'b/d':[0.1,0.6],
            'Lsh/d':[5,8],
            'b/Lsh':[0.002,0.05],
            'Lsh/b':[0.0,100],
            'Limp/d':[5,8]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        a6=0.44
        v_bl=a6*b*np.sqrt(Su*d/m*(42.7+Lsh/b))
        super().check('v_bl',v_bl)
        super().check('b/d',b/d)
        super().check('Lsh/d',Lsh/d)
        super().check('b/Lsh',b/Lsh)
        super().check('Lsh/b',Lsh/b)
        super().check('Limp/d',Limp/d)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
            b,d,m,v,Su,Lsh=symbols('b d m v Su Lsh')
            a6=0.4
            var('a6')
            g=a6*b*sqrt(Su*d/m*(42.7+Lsh/b))-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,v))
            self.d4=str(diff(g,Su))
            self.d5=str(diff(g,Lsh))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            Su=X[4]
            Lsh=X[5]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            Su=X[4]
            Lsh=X[5]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            dGdX[5] =eval(self.d5)
            super().SetdGdX(dGdX)
################################
#             SwRI              #
################################
class SwRI(penMed):
    """
    Southwest Research Institute (SwRI) model (Baker et al., 1980)
    ---
    b   thickness of a shield
    m   initial mass of the impactor
    v   velocity of impactor
    th  angle between a normal vector to a shield surface and the direction of impactor
    fragment  :'Standard' or 'Alternative'
    """
    S=0
    b1=0
    b2=0
    b3=0
    def __init__(self):
        self.variable=['b','m','v','th']
        self.title='Southwest Research Institute (SwRI) model'
        super().SaveRange('Validation process is not defined.')
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global S,b1,b2,b3
        tab={"0":{"b1":1414,"b2":0.295,"b3":0.910},
        "1":{"b1":1936,"b2":0.096,"b3":1.310},
        "2":{"b1":2039,"b2":0.064,"b3":0.430}}
        m=data['m']['mean']
        b=data['b']['mean']
        v=data['v']['mean']
        th=data['th']['mean']
        if data['fragment'] == 'Standard':
            k=0.186
        else:
            k=0.34
        S=1.33*(m/k)**(2/3)
        z=b/np.sqrt(S)
        if z>0 and z<=0.46:
            a="0"
        if z>0.46 and z<=1.06:
            a="1"
        if z>1.06:
            a="2"
        b1=tab[a]["b1"]
        b2=tab[a]["b2"]
        b3=tab[a]["b3"]
        print('Validation process is not defined.')
    class G(ls.Lbase):
        def __init__(self,n):
            global S,b1,b2,b3
            self.n=n
            super().__init__(self.n)
            b,m,v,th=symbols('b m v th')
            var('S b1 b2 b3')
            g=0.205*b1/sqrt(m)*S**b2*(39.37*b/cos(th))**b3-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,m))
            self.d2=str(diff(g,v))
            self.d4=str(diff(g,th))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            m=X[1]
            v=X[2]
            th=X[3]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            m=X[1]
            v=X[2]
            th=X[3]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            super().SetdGdX(dGdX)
################################
#             Lambert          #
################################
class Lambert(penMed):
    """
    Lambert and Jonas Approximation (1976)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    th  angle between a normal vector to a shield surface and the direction of impact
    v   velocity of impactor
    Limp    length of impactor
    ***constants***   
    ro_imp  material density of impactor
    a10     1750 for aluminum and 4000 for rolled homogeneous armor(RHA)
    *** 注意 ***
    SetMaterial(mat)で材料指定のこと
    mat: 'aluminum' or 'RHA'
          RHA:rolled homogeneous armor
    """
    a10=0
    def __init__(self):
        self.variable=['b','d','m','th','v','Limp']
        self.title='Lambert and Jonas Approximation'
        val_range={
            'm':[0.0005,3.63],
            'd':[0.002,0.05],
            'Limp/d':[4,30],
            'b':[0.006,0.15],
            'th':[0,60/360],
            'ro_imp':[7800,19000]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global a10
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        th=data['th']['mean']
        Limp=data['Limp']['mean']
        ro_imp=data['ro_imp']['mean']
        super().check('b',b)
        super().check('d',d)
        super().check('m',m)
        super().check('th',th)
        super().check('Limp/d',Limp/d)
        super().check('ro_imp',ro_imp)
    def SetMaterial(self,mat):
        global a10
        if mat=='aluminum':
            a10=1750
        if mat=='RHA':
            a10=4000
    class G(ls.Lbase):
        def __init__(self,n):
            global a10
            self.n=n
            super().__init__(self.n)
            b,d,m,th,v,Limp=symbols('b d m th v Limp')
            z=(b/d)*(1/cos(th))**0.75
            f=z+exp(z)-1
            var('a10')
            g=31.62*a10*(Limp/d)**0.15*sqrt(f*d**3/m)-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,th))
            self.d4=str(diff(g,v))
            self.d5=str(diff(g,Limp))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            th=X[3]
            v=X[4]
            Limp=X[5]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            th=X[3]
            v=X[4]
            Limp=X[5]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            dGdX[5] =eval(self.d5)
            super().SetdGdX(dGdX)
################################
#             Neilson          #
################################
class Neilson(penMed):
    """
    Neilson Formula (1985)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Su  ultimate tensile strength of shield material
    Lsh unsupported shield panel span
    v   velocity of impactor
    ***constants***
    Limp   length of impactor
    *** 注意 ***
    SetNose(Lsh,d,n_shape)で寸法、ノーズ形状指定のこと
    n_shape: 'flat' or 'hemispherical'
    """
    def __init__(self):
        self.variable=['b','d','m','Su','Lsh','v']
        self.title='Nelson Formula'
        val_range={
            'b/d':[0.14,0.64],
            'Lsh/d':[4,22],
            'Limp/d':[13,1e4],
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global a10,Limp
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        super().check('b/d',b/d)
        #super().check('Lsh/d',Lsh/d)
        super().check('Limp/d',Limp/d)
    def SetNose(self,Lsh,d,n_shape):
        global a12
        if n_shape=='flat':
            if Lsh/d > 4.0 and Lsh/d<22.0:
                a12=1.67
            if Lsh/d >=22.0:
                a12=4.26
        if n_shape=='hemispherical':
            a12=4.24
    class G(ls.Lbase):
        def __init__(self,n):
            global a12
            global a10,Limp
            self.n=n
            super().__init__(self.n)
            b,d,m,Su,Lsh,v=symbols('b d m Su Lsh v')
            #a12=1.67
            var('a12')
            g=a12*d*sqrt(Su*d/m)*(b/d)**0.85*(Lsh/d)**0.3-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,Su))
            self.d4=str(diff(g,Lsh))
            self.d5=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            dGdX[5] =eval(self.d5)
            super().SetdGdX(dGdX)
################################
#             Jowett         #
################################
class Jowett(penMed):
    """
    Jowett Formula (1986)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Su  ultimate tensile strength of shield material
    v   velocity of impactor
    ***constants***
    Lsh unsupported shield panel span
    Limp   length of impactor
    ***method***
    SetParam(b,d,Lsh)  g計算のためのparameter設定
    """
    def __init__(self):
        self.variable=['b','d','m','Su','v']
        self.title='Jowett Formula'
        val_range={
            'vbl':[40,200],
            'Su':[315,483],
            'Limp/d':[2,8],
            'b/d':[0.1,0.64]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def SetParam(self,b,d,Lsh):
        global ratio,omg
        if Lsh/d<=12:
            omg=(Lsh/d)**0.305
        else:
            #omg=12.0
            omg=2.14
        ratio=b/d       
    def Validation(self,data):
        global ratio,omg
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        if Lsh/d<=12:
            omg=(Lsh/d)**0.305
        else:
            #omg=12.0
            omg=2.14
        vbl=0
        if b/d>0.1 and b/d <0.25:
            vbl=1.62*omg*d*np.sqrt(Su*d/m)*(b/d)**0.87
        if b/d>=0.25 and b/d<0.64:
            vbl=0.87*omg*d*np.sqrt(Su*d/m)*(b/d)**0.42
        ratio=b/d
        super().check('vbl',vbl)
        super().check('Su',Su)
        super().check('Limp/d',Limp/d)
        super().check('b/d',b/d)
    class G(ls.Lbase):
        def __init__(self,n):
            global ratio,omg
            self.n=n
            super().__init__(self.n)
            b,d,m,Su,v=symbols('b d m Su v')
            vbl=0
            var('ratio omg')
            if ratio>0.1 and ratio <0.25:
                vbl=1.62*omg*d*sqrt(Su*d/m)*(b/d)**0.87
            if ratio>=0.25 and ratio<0.64:
                vbl=0.87*omg*d*sqrt(Su*d/m)*(b/d)**0.42
            g=vbl-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,Su))
            self.d4=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            v=X[4]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            v=X[4]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            super().SetdGdX(dGdX)
################################
#             WenJones         #
################################
class WenJones(penMed):
    """
    Wen and Jones Formula (1992)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Sy  yield stress of shield material
    Lsh unsupported shield panel span
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','Sy','Lsh','v']
        self.title='Jowett Formula'
        val_range={
            'vbl':[0,20],
            'Su':[340,440],
            'Lsh/d':[40,40],
            'Lsh/b':[25,100],
            'b/d':[0.4,1.6]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Sy=data['Sy']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        vbl=2*d*np.sqrt(Sy*d/m*(0.25*np.pi*(b/d)**2+(b/d)**1.47*(Lsh/d)**0.21))
        super().check('vbl',vbl)
        super().check('Su',Su)
        super().check('Lsh/d',Lsh/d)
        super().check('Lsh/b',Lsh/b)
        super().check('b/d',b/d)
    class G(ls.Lbase):
        def __init__(self,n):
            global ratio,omg
            self.n=n
            super().__init__(self.n)
            b,d,m,Sy,Lsh,v=symbols('b d m Sy Lsh v')
            g=2*d*sqrt(Sy*d/m*(0.25*np.pi*(b/d)**2+(b/d)**1.47*(Lsh/d)**0.21))-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,Sy))
            self.d4=str(diff(g,Lsh))
            self.d5=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Sy=X[3]
            Lsh=X[4]
            v=X[5]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Sy=X[3]
            Lsh=X[4]
            v=X[5]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            dGdX[5] =eval(self.d5)
            super().SetdGdX(dGdX)
################################
#             AlyLi         #
################################
class AlyLi(penMed):
    """
    Aly and Li Formulas (2008)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Su  ultimate tensile strength of shield material
    Lsh unsupported shield panel span
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','Su','Lsh','v']
        self.title='Aly and Li Formulas'
        val_range={
            'b/d':[0.1,0.64]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global v_Lsh,v_d,v_b
        v_b=data['b']['mean']
        v_d=data['d']['mean']
        v_Lsh=data['Lsh']['mean']
        super().check('b/d',v_b/v_d)
    class G(ls.Lbase):
        def __init__(self,n):
            global v_Lsh,v_d,v_b
            self.n=n
            super().__init__(self.n)
            b,d,m,Su,Lsh,v=symbols('b d m Su Lsh v')
            vbl=0
            var('v_Lsh v_d v_b')
            if v_Lsh/v_d <=12:
                if v_b/v_d > 0.1 and v_b/v_d<0.25:
                    vbl=1.79*d*sqrt(Su*d/m)*(b/d)**0.87*(Lsh/d)**0.305
                if v_b/v_d>=0.25 and v_b/v_d<0.64:
                    vbl=1.72*d*sqrt(Su*d/m)*(b/d)**0.42*(Lsh/d)**0.35
            if v_Lsh/v_d>12:
                if v_b/v_d > 0.1 and v_b/v_d<0.25:
                    vbl=3.44*d*sqrt(Su*d/m)*(b/d)**0.78
                if v_b/v_d>=0.25 and v_b/v_d<0.64:
                    vbl=1.72*d*sqrt(Su*d/m)*(b/d)**0.41                
            g=vbl-v
            self.gg=str(g)
            self.d0=str(diff(g,b))
            self.d1=str(diff(g,d))
            self.d2=str(diff(g,m))
            self.d3=str(diff(g,Su))
            self.d4=str(diff(g,Lsh))
            self.d5=str(diff(g,v))
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            g=eval(self.gg)
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            dGdX[0]=eval(self.d0)
            dGdX[1] =eval(self.d1)
            dGdX[2] =eval(self.d2)
            dGdX[3] =eval(self.d3)
            dGdX[4] =eval(self.d4)
            dGdX[5] =eval(self.d5)
            super().SetdGdX(dGdX)
################################
#             AlyLi_M         #
################################
class AlyLi_M(penMed):
    """
    Aly and Li Formulas (2008)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Su  ultimate tensile strength of shield material
    Lsh unsupported shield panel span
    v   velocity of impactor
    ---
    注意:解析前に必ずValidation(data)を実行すること
    """
    def __init__(self):
        self.variable=['b','d','m','Su','Lsh','v']
        self.title='Aly and Li Formulas'
        val_range={
            'b/d':[0.1,0.64]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global v_Lsh,v_d,v_b
        v_b=data['b']['mean']
        v_d=data['d']['mean']
        v_Lsh=data['Lsh']['mean']
        super().check('b/d',v_b/v_d)
    class G(ls.Lbase):
        def __init__(self,n):
            global v_Lsh,v_d,v_b
            self.n=n
            super().__init__(self.n)

        def gcalc(self):
            global v_Lsh,v_d,v_b
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            if v_Lsh/v_d <=12:
                if v_b/v_d > 0.1 and v_b/v_d<0.25:
                    vbl=1.79*d*np.sqrt(Su*d/m)*(b/d)**0.87*(Lsh/d)**0.305
                if v_b/v_d>=0.25 and v_b/v_d<0.64:
                    vbl=1.72*d*np.sqrt(Su*d/m)*(b/d)**0.42*(Lsh/d)**0.35
            if v_Lsh/v_d>12:
                if v_b/v_d > 0.1 and v_b/v_d<0.25:
                    vbl=3.44*d*np.sqrt(Su*d/m)*(b/d)**0.78
                if v_b/v_d>=0.25 and v_b/v_d<0.64:
                    vbl=1.72*d*np.sqrt(Su*d/m)*(b/d)**0.41                
            g=vbl-v
            super().SetG(g)
        def dGdXcalc(self):
            global v_Lsh,v_d,v_b
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            if v_Lsh/v_d <=12:
                if v_b/v_d > 0.1 and v_b/v_d<0.25:
                    dGdX[0]= 1.5573*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/b
                    dGdX[1]= 0.58175*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)
                    dGdX[2]= -0.895*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/m
                    dGdX[3]= 0.895*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/Su
                    dGdX[4]= 0.54595*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/Lsh
                    dGdX[5]= -1
                if v_b/v_d>=0.25 and v_b/v_d<0.64:
                    dGdX[0]= 0.7224*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/b
                    dGdX[1]= 1.2556*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)
                    dGdX[2]= -0.86*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/m
                    dGdX[3]= 0.86*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/Su
                    dGdX[4]= 0.602*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/Lsh
                    dGdX[5]= -1
            if v_Lsh/v_d>12:
                if v_b/v_d > 0.1 and v_b/v_d<0.25:
                    dGdX[0]= 2.6832*d*(b/d)**0.78*np.qrt(Su*d/m)/b
                    dGdX[1]= 2.4768*(b/d)**0.78*np.sqrt(Su*d/m)
                    dGdX[2]= -1.72*d*(b/d)**0.78*np.sqrt(Su*d/m)/m
                    dGdX[3]= 1.72*d*(b/d)**0.78*np.sqrt(Su*d/m)/Su
                    dGdX[4]= 0
                    dGdX[5]= -1
                if v_b/v_d>=0.25 and v_b/v_d<0.64:
                    dGdX[0]= 0.7052*d*(b/d)**0.41*np.sqrt(Su*d/m)/b
                    dGdX[1]= 1.8748*(b/d)**0.41*np.sqrt(Su*d/m)
                    dGdX[2]= -0.86*d*(b/d)**0.41*np.sqrt(Su*d/m)/m
                    dGdX[3]= 0.86*d*(b/d)**0.41*np.sqrt(Su*d/m)/Su
                    dGdX[4]= 0
                    dGdX[5]= -1  
            super().SetdGdX(dGdX)
################################
#             THOR_M             #
################################
class THOR_M(penMed):
    """
    THOR equation (Crull and Swisdak, 2005)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    th  angle between a normal vector to a shield surface and the direction of impactor
    """
    C1=0
    a1=0
    b1=0
    g1=0
    def __init__(self):
        self.variable=['b','d','m','v','th']
        self.title='THOR Equations'
        self.tab={"Magnesium":{"C1":6.349,"a1":1.004,"b1":-1.076,"g1":0.966},
             "Aluminum":{"C1":6.185,"a1":0.903,"b1":-0.941,"g1":1.098},
             "CastIron":{"C1":10.153,"a1":2.186,"b1":-2.204,"g1":2.156},
             "Titanium":{"C1":7.552,"a1":1.325,"b1":-1.314,"g1":1.643},
             "FaceSteel":{"C1":7.694,"a1":1.191,"b1":-1.397,"g1":1.747}, 
             "MildSteel":{"C1":6.523,"a1":0.906,"b1":-0.963,"g1":1.286}, 
             "HardSteel":{"C1":6.601,"a1":0.906,"b1":-0.963,"g1":1.286},
             "Copper":{"C1":14.065,"a1":3.476,"b1":-3.687,"g1":4.27},
             "Lead":{"C1":10.955,"a1":2.735,"b1":-2.753,"g1":3.59}
            }
        super().SaveRange('Validation process is not defined.')
        super().SaveVariable(self.variable)
    def MatList(self):
        """
        目的:登録されている材料名リストを返す
        """
        return list(self.tab.keys())
    def setMaterial(self,mat):
        """
        目的:材料の設定
           mat: "Magnesium","Aluminum","CastIron","Titanium","FaceSteel","MildSteel","HardSteel","Copper","Lead"
        """
        global C1,a1,b1,g1
        C1=self.tab[mat]["C1"]
        a1=self.tab[mat]["a1"]
        b1=self.tab[mat]["b1"]
        g1=self.tab[mat]["g1"]
    class G(ls.Lbase):
        def __init__(self,n):
            global C1,a1,b1,g1
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            global C1,a1,b1,g1
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            th=X[4]
            A=np.pi*d*d/4
            g=0.3048*10**C1*(61024*b*A)**a1*(15432.4*m)**b1*(1/np.cos(th))**g1-v
            super().SetG(g)
        def dGdXcalc(self):
            global C1,a1,b1,g1
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            th=X[4]
            dGdX[0]=0.3048*10**C1*a1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1/b
            dGdX[1] =0.6096*10**C1*a1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1/d
            dGdX[2] =0.3048*10**C1*b1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1/m
            dGdX[3] =-1
            dGdX[4] =0.3048*10**C1*g1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1*np.sin(th)/np.cos(th)
            super().SetdGdX(dGdX)
################################
#             BRL_M            #
################################
class BRL_M(penMed):
    """
    Ballistic Research Laboratories (BRL) model (1968)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    ***constants***
    Limp    length of impactor
    Lsh     unsupported shield panel span
    Su      ultimate tensile strength of shield material
    """
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='BRL Formula'
        val_range={
            'v_bl':[57,270],
            'Limp/d':[1.25,8],
            'b/d':[0.1,1.0],
            'Lsh/d':[8,35],
            'Su':[315,500]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Limp=data['Limp']['mean']
        Lsh=data['Lsh']['mean']
        Su=data['Su']['mean']
        a7=5.37
        v_bl=a7*1e4*(b*d)**0.75/m**0.5
        super().check('v_bl',v_bl)
        super().check('Limp/d',Limp/d)
        super().check('b/d',b/d)
        super().check('Lsh/d',Lsh/d)
        super().check('Su',Su)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            a7=5.37
            g=a7*1e4*(b*d)**0.75/m**0.5-v
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            a7=5.37
            dGdX[0]= 7500.0*a7*(b*d)**0.75/(b*m**0.5)
            dGdX[1]= 7500.0*a7*(b*d)**0.75/(d*m**0.5)
            dGdX[2]= -5000.0*a7*(b*d)**0.75/m**1.5
            dGdX[3] =-1
            super().SetdGdX(dGdX)
################################
#             DeMarre_M        #
################################
class DeMarre_M(penMed):
    """
    De Marre formula (Herrmann and Jones,1961)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='De Marre Formula'
        val_range={
            'v_bl':[200,900],
            'm':[0.1,50]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        v_bl=0.4311e5*d**0.75*b**0.7/m**0.5
        super().check('v_bl',v_bl)
        super().check('m',m)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)

        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            g=0.431e5*d**0.75*b**0.7/m**0.5-v
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            dGdX[0]=30170.0*d**0.75/(b**0.3*m**0.5)
            dGdX[1] =32325.0*b**0.7/(d**0.25*m**0.5)
            dGdX[2] =-21550.0*b**0.7*d**0.75/m**1.5
            dGdX[3] =-1
            super().SetdGdX(dGdX)
################################
#             Jowett_M         #
################################
class Jowett_M(penMed):
    """
    Jowett Formula (1986)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Su  ultimate tensile strength of shield material
    v   velocity of impactor
    ***constants***
    Lsh unsupported shield panel span
    Limp   length of impactor
    ***method***
    SetParam(b,d,Lsh)  g計算のためのparameter設定
    """
    def __init__(self):
        self.variable=['b','d','m','Su','v']
        self.title='Jowett Formula'
        val_range={
            'vbl':[40,200],
            'Su':[315,483],
            'Limp/d':[2,8],
            'b/d':[0.1,0.64]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def SetParam(self,b,d,Lsh):
        global ratio,omg
        if Lsh/d<=12:
            omg=(Lsh/d)**0.305
        else:
            #omg=12.0
            omg=2.14
        ratio=b/d       
    def Validation(self,data):
        global ratio,omg
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        if Lsh/d<=12:
            omg=(Lsh/d)**0.305
        else:
            #omg=12.0
            omg=2.14
        vbl=0
        if b/d>0.1 and b/d <0.25:
            vbl=1.62*omg*d*np.sqrt(Su*d/m)*(b/d)**0.87
        if b/d>=0.25 and b/d<0.64:
            vbl=0.87*omg*d*np.sqrt(Su*d/m)*(b/d)**0.42
        ratio=b/d
        super().check('vbl',vbl)
        super().check('Su',Su)
        super().check('Limp/d',Limp/d)
        super().check('b/d',b/d)
    class G(ls.Lbase):
        def __init__(self,n):
            global ratio,omg
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            global ratio,omg
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            v=X[4]
            ratio=b/d
            vbl=0
            if ratio>0.1 and ratio <0.25:
                vbl=1.62*omg*d*np.sqrt(Su*d/m)*(b/d)**0.87
            if ratio>=0.25 and ratio<0.64:
                vbl=0.87*omg*d*np.sqrt(Su*d/m)*(b/d)**0.42
            g=vbl-v
            super().SetG(g)
        def dGdXcalc(self):
            global ratio,omg
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            v=X[4]
            ratio=b/d
            if ratio>0.1 and ratio <0.25:
                dGdX[0]= 1.4094*d*omg*(b/d)**0.87*np.sqrt(Su*d/m)/b
                dGdX[1]= 1.0206*omg*(b/d)**0.87*np.sqrt(Su*d/m)
                dGdX[2]= -0.81*d*omg*(b/d)**0.87*np.sqrt(Su*d/m)/m
                dGdX[3]= 0.81*d*omg*(b/d)**0.87*np.sqrt(Su*d/m)/Su
                dGdX[4]= -1
            if ratio>=0.25 and ratio<0.64:
                dGdX[0]= 0.3654*d*omg*(b/d)**0.42*np.sqrt(Su*d/m)/b
                dGdX[1]= 0.9396*omg*(b/d)**0.42*np.sqrt(Su*d/m)
                dGdX[2]= -0.435*d*omg*(b/d)**0.42*np.sqrt(Su*d/m)/m
                dGdX[3]= 0.435*d*omg*(b/d)**0.42*np.sqrt(Su*d/m)/Su
                dGdX[4]= -1
            super().SetdGdX(dGdX)
################################
#             Lambert_M        #
################################
class Lambert_M(penMed):
    """
    Lambert and Jonas Approximation (1976)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    th  angle between a normal vector to a shield surface and the direction of impact
    v   velocity of impactor
    Limp    length of impactor
    ***constants***   
    ro_imp  material density of impactor
    a10     1750 for aluminum and 4000 for rolled homogeneous armor(RHA)
    *** 注意 ***
    SetMaterial(mat)で材料指定のこと
    mat: 'aluminum' or 'RHA'
          RHA:rolled homogeneous armor
    """
    a10=0
    def __init__(self):
        self.variable=['b','d','m','th','v','Limp']
        self.title='Lambert and Jonas Approximation'
        val_range={
            'm':[0.0005,3.63],
            'd':[0.002,0.05],
            'Limp/d':[4,30],
            'b':[0.006,0.15],
            'th':[0,60/360],
            'ro_imp':[7800,19000]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global a10
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        th=data['th']['mean']
        Limp=data['Limp']['mean']
        ro_imp=data['ro_imp']['mean']
        super().check('b',b)
        super().check('d',d)
        super().check('m',m)
        super().check('th',th)
        super().check('Limp/d',Limp/d)
        super().check('ro_imp',ro_imp)
    def SetMaterial(self,mat):
        global a10
        if mat=='aluminum':
            a10=1750
        if mat=='RHA':
            a10=4000
    class G(ls.Lbase):
        def __init__(self,n):
            global a10
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            global a10
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            th=X[3]
            v=X[4]
            Limp=X[5]
            z=(b/d)*(1/np.cos(th))**0.75
            f=z+np.exp(z)-1
            g=31.62*a10*(Limp/d)**0.15*np.sqrt(f*d**3/m)-v
            super().SetG(g)
        def dGdXcalc(self):
            global a10
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            th=X[3]
            v=X[4]
            Limp=X[5]
            dGdX[0]= 15.81*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/cos(th))**0.75/d) - 1)/m)*((1/np.cos(th))**0.75*np.exp(b*(1/np.cos(th))**0.75/d)/d + (1/np.cos(th))**0.75/d)/(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)
            dGdX[1]= -4.743*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)/d + 31.62*a10*m*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)*(d**3*(-b*(1/np.cos(th))**0.75*np.exp(b*(1/np.cos(th))**0.75/d)/d**2 - b*(1/np.cos(th))**0.75/d**2)/(2*m) + 3*d**2*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/(2*m))/(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1))
            dGdX[2]= -15.81*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)/m
            dGdX[3]= 15.81*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)*(0.75*b*(1/np.cos(th))**0.75*np.exp(b*(1/np.cos(th))**0.75/d)*np.sin(th)/(d*np.cos(th)) + 0.75*b*(1/np.cos(th))**0.75*np.sin(th)/(d*np.cos(th)))/(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)
            dGdX[4]= -1
            dGdX[5]= 4.743*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)/Limp
            super().SetdGdX(dGdX)
################################
#             Neilson_M        #
################################
class Neilson_M(penMed):
    """
    Neilson Formula (1985)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Su  ultimate tensile strength of shield material
    Lsh unsupported shield panel span
    v   velocity of impactor
    ***constants***
    Limp   length of impactor
    *** 注意 ***
    SetNose(Lsh,d,n_shape)で寸法、ノーズ形状指定のこと
    n_shape: 'flat' or 'hemispherical'
    """
    def __init__(self):
        self.variable=['b','d','m','Su','Lsh','v']
        self.title='Nelson Formula'
        val_range={
            'b/d':[0.14,0.64],
            'Lsh/d':[4,22],
            'Limp/d':[13,1e4],
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global a10,Limp
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        super().check('b/d',b/d)
        #super().check('Lsh/d',Lsh/d)
        super().check('Limp/d',Limp/d)
    def SetNose(self,Lsh,d,n_shape):
        global a12
        if n_shape=='flat':
            if Lsh/d > 4.0 and Lsh/d<22.0:
                a12=1.67
            if Lsh/d >=22.0:
                a12=4.26
        a12=1.67 #上の条件を無視
        if n_shape=='hemispherical':
            a12=4.24
    class G(ls.Lbase):
        def __init__(self,n):
            global a12
            global a10,Limp
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            global a12
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            g=a12*d*np.sqrt(Su*d/m)*(b/d)**0.85*(Lsh/d)**0.3-v
            super().SetG(g)
        def dGdXcalc(self):
            global a12
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Su=X[3]
            Lsh=X[4]
            v=X[5]
            dGdX[0]= 0.85*a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/b
            dGdX[1]= 0.35*a12*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)
            dGdX[2]= -a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/(2*m)
            dGdX[3]= a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/(2*Su)
            dGdX[4]= 0.3*a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/Lsh
            dGdX[5]= -1
            super().SetdGdX(dGdX)
################################
#             Ohte_M           #
################################
class Ohte_M(penMed):
    """
    Ohte et al. Formula (Ohte et al., 1982)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='Ohte et al. Formula'
        val_range={
            'v_bl':[25,180],
            'm':[3,50],
            'Su':[490,637],
            'b':[7,38],
            'Lsh/b':[39,1e4],
            'd':[39,1e4]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        v_bl=7.67e4*(b*d)**0.75/m**0.5
        super().check('v_bl',v_bl)
        super().check('m',m)
        super().check('Su',Su)
        super().check('b',b)
        super().check('Lsh/b',Lsh/b)
        super().check('d',d)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            g=7.67e4*(b*d)**0.75/m**0.5-v
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            dGdX[0]= 57525.0*(b*d)**0.75/(b*m**0.5)
            dGdX[1]= 57525.0*(b*d)**0.75/(d*m**0.5)
            dGdX[2]= -38350.0*(b*d)**0.75/m**1.5
            dGdX[5]= -1
            super().SetdGdX(dGdX)
################################
#             SRI_M            #
################################
class SRI_M(penMed):
    """
    Stanford Research Institute (SRI) correlation (1963)
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    Lsh unsupported shield panel span
    Su  ultimate tensile strength of shield material
    Limp    length of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','v','Lsh','Su']
        self.title='Ohte et al. Formula'
        val_range={
            'v_bl':[21,122],
            'b/d':[0.1,0.6],
            'Lsh/d':[5,8],
            'b/Lsh':[0.002,0.05],
            'Lsh/b':[0.0,100],
            'Limp/d':[5,8]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        a6=0.44
        v_bl=a6*b*np.sqrt(Su*d/m*(42.7+Lsh/b))
        super().check('v_bl',v_bl)
        super().check('b/d',b/d)
        super().check('Lsh/d',Lsh/d)
        super().check('b/Lsh',b/Lsh)
        super().check('Lsh/b',Lsh/b)
        super().check('Limp/d',Limp/d)
    class G(ls.Lbase):
        def __init__(self,n):
            global a6
            self.n=n
            super().__init__(self.n)
            a6=0.4
        def gcalc(self):
            global a6
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            Su=X[4]
            Lsh=X[5]
            g=a6*b*np.sqrt(Su*d/m*(42.7+Lsh/b))-v
            super().SetG(g)
        def dGdXcalc(self):
            global a6
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            v=X[3]
            Su=X[4]
            Lsh=X[5]
            dGdX[0]= -Lsh*a6*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*b*(Lsh/b + 42.7)) + a6*np.sqrt(Su*d*(Lsh/b + 42.7)/m)
            dGdX[1]= a6*b*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*d)
            dGdX[2]= -a6*b*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*m)
            dGdX[3]= -1
            dGdX[4]= a6*b*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*Su)
            dGdX[5]= a6*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*(Lsh/b + 42.7))
            super().SetdGdX(dGdX)
################################
#             SwRI_M           #
################################
class SwRI_M(penMed):
    """
    Southwest Research Institute (SwRI) model (Baker et al., 1980)
    ---
    b   thickness of a shield
    m   initial mass of the impactor
    v   velocity of impactor
    th  angle between a normal vector to a shield surface and the direction of impactor
    fragment  :'Standard' or 'Alternative'
    """
    S=0
    b1=0
    b2=0
    b3=0
    def __init__(self):
        self.variable=['b','m','v','th']
        self.title='Southwest Research Institute (SwRI) model'
        super().SaveRange('Validation process is not defined.')
        super().SaveVariable(self.variable)
    def Validation(self,data):
        global S,b1,b2,b3
        tab={"0":{"b1":1414,"b2":0.295,"b3":0.910},
        "1":{"b1":1936,"b2":0.096,"b3":1.310},
        "2":{"b1":2039,"b2":0.064,"b3":0.430}}
        m=data['m']['mean']
        b=data['b']['mean']
        v=data['v']['mean']
        th=data['th']['mean']
        if data['fragment'] == 'Standard':
            k=0.186
        else:
            k=0.34
        S=1.33*(m/k)**(2/3)
        z=b/np.sqrt(S)
        if z>0 and z<=0.46:
            a="0"
        if z>0.46 and z<=1.06:
            a="1"
        if z>1.06:
            a="2"
        b1=tab[a]["b1"]
        b2=tab[a]["b2"]
        b3=tab[a]["b3"]
        print('Validation process is not defined.')
    class G(ls.Lbase):
        def __init__(self,n):
            global S,b1,b2,b3
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            global S,b1,b2,b3
            X=super().GetX()
            b=X[0]
            m=X[1]
            v=X[2]
            th=X[3]
            g=0.205*b1/np.sqrt(m)*S**b2*(39.37*b/cos(th))**b3-v
            super().SetG(g)
        def dGdXcalc(self):
            global S,b1,b2,b3
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            m=X[1]
            v=X[2]
            th=X[3]
            dGdX[0]= 0.205*S**b2*b1*b3*(39.37*b/np.cos(th))**b3/(b*np.sqrt(m))
            dGdX[1]= -0.1025*S**b2*b1*(39.37*b/np.cos(th))**b3/m**(3/2)
            dGdX[2]= -1
            dGdX[3]= 0.205*S**b2*b1*b3*(39.37*b/np.cos(th))**b3*np.sin(th)/(sqrt(m)*np.cos(th))
            super().SetdGdX(dGdX)
################################
#             WenJones_M       #
################################
class WenJones_M(penMed):
    """
    Wen and Jones Formula (1992)
    ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    Sy  yield stress of shield material
    Lsh unsupported shield panel span
    v   velocity of impactor
    """
    def __init__(self):
        self.variable=['b','d','m','Sy','Lsh','v']
        self.title='Jowett Formula'
        val_range={
            'vbl':[0,20],
            'Su':[340,440],
            'Lsh/d':[40,40],
            'Lsh/b':[25,100],
            'b/d':[0.4,1.6]
        }
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Sy=data['Sy']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        vbl=2*d*np.sqrt(Sy*d/m*(0.25*np.pi*(b/d)**2+(b/d)**1.47*(Lsh/d)**0.21))
        super().check('vbl',vbl)
        super().check('Su',Su)
        super().check('Lsh/d',Lsh/d)
        super().check('Lsh/b',Lsh/b)
        super().check('b/d',b/d)
    class G(ls.Lbase):
        def __init__(self,n):
            global ratio,omg
            self.n=n
            super().__init__(self.n)
        def gcalc(self):
            X=super().GetX()
            b=X[0]
            d=X[1]
            m=X[2]
            Sy=X[3]
            Lsh=X[4]
            v=X[5]
            g=2*d*np.sqrt(Sy*d/m*(0.25*np.pi*(b/d)**2+(b/d)**1.47*(Lsh/d)**0.21))-v
            super().SetG(g)
        def dGdXcalc(self):
            X=super().GetX()
            dGdX=super().GetdGdX()
            b=X[0]
            d=X[1]
            m=X[2]
            Sy=X[3]
            Lsh=X[4]
            v=X[5]
            dGdX[0]= d*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)*(1.5707963267949*b/d**2 + 1.47*(Lsh/d)**0.21*(b/d)**1.47/b)/(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)
            dGdX[1]= 2*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m) + 2*m*sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)*(Sy*d*(-1.5707963267949*b**2/d**3 - 1.68*(Lsh/d)**0.21*(b/d)**1.47/d)/(2*m) + Sy*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/(2*m))/(Sy*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47))
            dGdX[2]= -d*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)/m
            dGdX[3]= d*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)/Sy
            dGdX[4]= 0.21*d*(Lsh/d)**0.21*(b/d)**1.47*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)/(Lsh*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47))
            dGdX[5]= -1
            super().SetdGdX(dGdX)