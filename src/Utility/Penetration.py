from Utility import LimitState as ls
from sympy import *
class penMed(ls.RelBase):
    """
    目的:貫通評価クラスが継承する直前のクラス
    """
    def PoFcontour(self,lmsf,data,x,y,xlabel,ylabel):
        """
        目的:JSONデータdata内の二つのパラメータxlabel,ylabelに関するPOF等高線データ作成
            lmsf:   貫通評価クラスのインスタンス
            data:   貫通評価計算のための入力データ(JSON形式)
            x:      x = np.arange(0.006, 0.02, 0.0002)などで発生するデータ列
            y:      y = np.arange(0.006, 0.02, 0.0002)などで発生するデータ列
            xlabel: data内の対象ラベル
            ylabel: data内の対象ラベル
        戻り値: X,Y,Z
            plt.pcolormesh(X, Y, Z, cmap='hsv')などにより等高線描画する
        """
        X, Y = np.meshgrid(x, y)
        ZZ=[]
        for iy in range(len(y)):    
            yy=Y[:,0][iy]
            data[ylabel]['mean']=yy
            za=[]
            for ix in range(len(x)):
                xx=X[0][ix]
                data[xlabel]['mean']=xx
                lmsf.Reliability(data)
                pof=lmsf.GetPOF()
                za.append(pof)
            ZZ.append(za)
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
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
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
        super().SaveRange('Validation process is not defined.')
    def setMaterial(self,mat):
        """
        目的:材料の設定
           mat: "Magnesium","Aluminum","CastIron","Titanium","FaceSteel","MildSteel","HardSteel","Copper","Lead"
        """
        global C1,a1,b1,g1
        tab={"Magnesium":{"C1":6.349,"a1":1.004,"b1":-1.076,"g1":0.966},
             "Aluminum":{"C1":6.185,"a1":0.903,"b1":-0.941,"g1":1.098},
             "CastIron":{"C1":10.153,"a1":2.186,"b1":-2.204,"g1":2.156},
             "Titanium":{"C1":7.552,"a1":1.325,"b1":-1.314,"g1":1.643},
             "FaceSteel":{"C1":7.694,"a1":1.191,"b1":-1.397,"g1":1.747}, 
             "MildSteel":{"C1":6.523,"a1":0.906,"b1":-0.963,"g1":1.286}, 
             "HardSteel":{"C1":6.601,"a1":0.906,"b1":-0.963,"g1":1.286},
             "Copper":{"C1":14.065,"a1":3.476,"b1":-3.687,"g1":4.27},
             "Lead":{"C1":10.955,"a1":2.735,"b1":-2.753,"g1":3.59}
            }
        C1=tab[mat]["C1"]
        a1=tab[mat]["a1"]
        b1=tab[mat]["b1"]
        g1=tab[mat]["g1"]
    class G(ls.Lbase):
        def __init__(self,n):
            global C1,a1,b1,g1
            self.n=n
            super().__init__(self.n)
            b,d,m,v,th=symbols('b d m v th')
            A=np.pi*d*d/4
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
            g=0.205*b1/np.sqrt(m)*S**b2*(39.37*b/cos(th))**b3-v
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