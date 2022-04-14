from Utility import LimitState as ls
from Utility import UnitConv as uc
from sympy import *
################################
#             BRL              #
################################
class BRL(ls.RelBase,uc.Rcheck):
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='BRL Formula'
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Limp=data['Limp']['mean']
        Lsh=data['Lsh']['mean']
        Su=data['Su']['mean']
        a7=5.37
        v_bl=a7*1e4*(b*d)**0.75/m**0.5
        super().check('v_bl',v_bl,57,270)
        super().check('Limp/d',Limp/d,1.25,8)
        super().check('b/d',b/d,0.1,1.0)
        super().check('Lsh/d',Lsh/d,8,35)
        super().check('Su',Su,315,500)
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
class DeMarre(ls.RelBase,uc.Rcheck):
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='De Marre Formula'
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        v_bl=0.4311e5*d**0.75*b**0.7/m**0.5
        super().check('v_bl',v_bl,200,900)
        super().check('m',m,0.1,50)
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

class THOR(ls.RelBase,uc.Rcheck):
    C1=0
    a1=0
    b1=0
    g1=0
    def __init__(self):
        self.variable=['b','d','m','v','th']
        self.title='THOR Equations'
    def setMaterial(self,mat):
        """
        目的:材料の設定
           mat: "Magnesium","Aluminum","CastIron","Titanium","FaceSteel","MildSteel","HardSteel","Copper","Lead"
        """
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
class Ohte(ls.RelBase,uc.Rcheck):
    def __init__(self):
        self.variable=['b','d','m','v']
        self.title='Ohte et al. Formula'
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        v_bl=7.67e4*(b*d)**0.75/m**0.5
        super().check('v_bl',v_bl,25,180)
        super().check('m',m,3,50)
        super().check('Su',Su,490,637)
        super().check('b',b,7,38)
        super().check('Lsh/b',Lsh/b,39,1e4)
        super().check('d',d,39,1e4)
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
class SRI(ls.RelBase,uc.Rcheck):
    def __init__(self):
        self.variable=['b','d','m','v','Lsh','Su']
        self.title='Ohte et al. Formula'
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        a6=0.44
        v_bl=a6*b*np.sqrt(Su*d/m*(42.7+Lsh/b))
        super().check('v_bl',v_bl,21,122)
        super().check('b/d',b/d,0.1,0.6)
        super().check('Lsh/d',Lsh/d,5,8)
        super().check('b/Lsh',b/Lsh,0.002,0.05)
        super().check('Lsh/b',Lsh/b,0.0,100)
        super().check('Limp/d',Limp/d,5,8)
    class G(ls.Lbase):
        def __init__(self,n):
            self.n=n
            super().__init__(self.n)
            b,d,m,v,Su,Lsh=symbols('b d m v Su Lsh')
            a6=0.44
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