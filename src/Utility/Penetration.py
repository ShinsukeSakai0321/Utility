from LimitState import LimitState as ls
import numpy as np
import pickle
class penMed(ls.LSFM):
    """
    目的:貫通評価クラスが継承する直前のクラス
    """
    def __init__(self):
        self.range={}
        self.variable=[]
        self.const=[]
        self.title=''
        self.df={}
    def MakeContour(self,data,cdata):
        """
        目的:JSONデータdata内の二つのパラメータに関するPOF等高線データ作成
            data    貫通評価計算のための入力データ(JSON形式)を基本とし、等高線を作成する
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
        ZZ=[[] for i in range(len(self.variable)+2)]
        for iy in range(len(y)):    
            yy=Y[:,0][iy]
            data[key[1]]['mean']=yy
            zPoF=[]
            zBeta=[]
            zAlpha=[[] for i in range(len(self.variable))]
            for ix in range(len(x)):
                xx=X[0][ix]
                data[key[0]]['mean']=xx
                self.df=data
                muX,sigmmaX,dist=self.makeData(data)
                super().__init__(len(muX),muX,sigmmaX,dist)
                super().RFn()
                zPoF.append(super().GetPOF())
                zBeta.append(super().GetBeta())
                for ii in range(len(self.variable)):
                    zAlpha[ii].append(super().GetAlpha()[ii])
            ZZ[0].append(zPoF)
            ZZ[1].append(zBeta)
            for ii in range(len(self.variable)):
                ZZ[ii+2].append(zAlpha[ii])
        return X,Y,ZZ
    def Calc(self,df):
        self.df=df
        muX,sigmmaX,dist=self.makeData(df)
        super().__init__(len(muX),muX,sigmmaX,dist)
        super().RFn()
    def SaveDf(self,df):
        self.df=df
    def gDict(self):
        """
        計算用の入力辞書データの呼び出し
        """
        return self.df
    def pklRead(self,fname):
        with open(fname,'rb') as f:
            df = pickle.load(f)
        return df
    def makeData(self,df):
        """
        辞書型データのLSFM入力データへの変換
        """
        muX=[]
        sigmmaX=[]
        dist=[]
        for var in self.variable:
            muX.append(df[var]['mean'])
            sigmmaX.append(df[var]['mean']*df[var]['cov'])
            dist.append(df[var]['dist'])
        return muX,sigmmaX,dist
    def makeMean(self,df):
        """
        辞書型データのmeanのみ抽出してもどす
        """
        muX=[]
        for var in self.variable:
            muX.append(df[var]['mean'])
        return muX                       
    def SaveRange(self,aa):
        """
        目的:適用範囲データの保存
        """
        self.Range=aa
    def SaveConst(self,aa):
        self.const=aa
    def SaveVariable(self,aa):
        self.variable=aa
    def ShowRange(self):
        """
        目的:適用範囲データの表示
        """
        return self.Range
    def SaveTitle(self,title):
        self.title=title
    def GetTitle(self):
        return self.title
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
    def check_c(self,cond,val):
        """
        目的:適用範囲内か否かについて、True,Falseでもどす
        """ 
        min_r=self.Range[cond][0]
        max_r=self.Range[cond][1]
        if val >= min_r and val<=max_r:
            iflag=True
        else:
            iflag=False
        return iflag       
    def GetMean(self,data):
        """
        目的:Jsonデータ　dataから、変数の平均値を取り出しリストを戻す
        """
        dmean=[]
        for aa in data.keys():
            if 'cov' in data[aa].keys() or 'sd' in data[aa].keys():
                dmean.append(data[aa]['mean'])
        return dmean
    def Validation(self,data):
        print('Varidation process is not defined')
################################
#             AlyLi         #
################################
#  以下_Mがつくクラスでは、数式処理を使わない
#
class AlyLi(penMed):
    """
    ---  Aly and Li Formulas (2008)  ---
    SY Aly and QM Li. ,Critical impact energy for the perforation
    of metallic plates. Nuclear
    Engineering and Design, Vol. 238, No. 10, pp. 2521–2528, 2008.
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
        self.variable=['b','d','m','Su','Lsh','v','Me']
        self.const=[]
        title='Aly and Li Formulas'
        val_range={
            'b/d':[0.1,0.64]
        }
        self.i_Valid=True  #Validation結果を出力するときTrue
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
        return
    def Validation(self,data):
        v_b=data['b']['mean']
        v_d=data['d']['mean']
        v_Lsh=data['Lsh']['mean']
        if self.i_Valid:
            super().check('b/d',v_b/v_d)
        else:
            ii=0
            if super().check_c('b/d',v_b/v_d)!=True: ii+=1
            return ii

    def gcalc(self,X):
        df=super().gDict()
        v_Lsh=df['Lsh']['mean']
        v_d=df['d']['mean']
        v_b=df['b']['mean']
        b=X[0]
        d=X[1]
        m=X[2]
        Su=X[3]
        Lsh=X[4]
        v=X[5]
        Me=X[6]
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
        g=Me*vbl-v
        return g
    def dGdXcalc(self,X):
        df=super().gDict()
        v_Lsh=df['Lsh']['mean']
        v_d=df['d']['mean']
        v_b=df['b']['mean']
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        Su=X[3]
        Lsh=X[4]
        v=X[5]
        Me=X[6]
        if v_Lsh/v_d <=12:
            if v_b/v_d > 0.1 and v_b/v_d<0.25:
                dGdX[0]= Me*1.5573*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/b
                dGdX[1]= Me*0.58175*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)
                dGdX[2]= Me*(-0.895)*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/m
                dGdX[3]= Me*0.895*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/Su
                dGdX[4]= Me*0.54595*d*(Lsh/d)**0.305*(b/d)**0.87*np.sqrt(Su*d/m)/Lsh
                dGdX[5]= -1
                dGdX[6]=1.79*d*np.sqrt(Su*d/m)*(b/d)**0.87*(Lsh/d)**0.305
            if v_b/v_d>=0.25 and v_b/v_d<0.64:
                dGdX[0]= Me*0.7224*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/b
                dGdX[1]= Me*1.2556*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)
                dGdX[2]= Me*(-0.86)*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/m
                dGdX[3]= Me*0.86*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/Su
                dGdX[4]= Me*0.602*d*(Lsh/d)**0.35*(b/d)**0.42*np.sqrt(Su*d/m)/Lsh
                dGdX[5]= -1
                dGdX[6]=1.72*d*np.sqrt(Su*d/m)*(b/d)**0.42*(Lsh/d)**0.35
        if v_Lsh/v_d>12:
            if v_b/v_d > 0.1 and v_b/v_d<0.25:
                dGdX[0]= Me*2.6832*d*(b/d)**0.78*np.sqrt(Su*d/m)/b
                dGdX[1]= Me*2.4768*(b/d)**0.78*np.sqrt(Su*d/m)
                dGdX[2]= Me*(-1.72)*d*(b/d)**0.78*np.sqrt(Su*d/m)/m
                dGdX[3]= Me*1.72*d*(b/d)**0.78*np.sqrt(Su*d/m)/Su
                dGdX[4]= 0
                dGdX[5]= -1
                dGdX[6]=3.44*d*np.sqrt(Su*d/m)*(b/d)**0.78
            if v_b/v_d>=0.25 and v_b/v_d<0.64:
                dGdX[0]= Me*0.7052*d*(b/d)**0.41*np.sqrt(Su*d/m)/b
                dGdX[1]= Me*1.8748*(b/d)**0.41*np.sqrt(Su*d/m)
                dGdX[2]= Me*(-0.86)*d*(b/d)**0.41*np.sqrt(Su*d/m)/m
                dGdX[3]= Me*0.86*d*(b/d)**0.41*np.sqrt(Su*d/m)/Su
                dGdX[4]= 0
                dGdX[5]= -1
                dGdX[6] =1.72*d*np.sqrt(Su*d/m)*(b/d)**0.41    
        return dGdX
################################
#             THOR             #
################################
class THOR(penMed):
    """
    ---  THOR equation (Crull and Swisdak, 2005)  ---
    M. Crull and Jr. Swisdak, M.M. Methodologies for calculating
    primary fragment charcteristics. Technical Paper No.16
    (Technical Report DDESB TP16),Revision 2. Department
    of Defense Explosives, Safety Board, Alexandria, VA, 2005.
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
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','v','th','Me']
        self.const=['Material']
        title='THOR Equations'
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
        super().SaveTitle(title)
        super().SaveRange('Validation process is not defined.')
        super().SaveVariable(self.variable)
    def Validation(self,data):
        if self.i_Valid:
            print('Validation process is not defined.')
        else:
            return 0
    def MatList(self):
        """
        目的:登録されている材料名リストを返す
        """
        return list(self.tab.keys())
    def gcalc(self,X):
        df=super().gDict()
        mat=df['Material']
        self.C1=self.tab[mat]["C1"]
        self.a1=self.tab[mat]["a1"]
        self.b1=self.tab[mat]["b1"]
        self.g1=self.tab[mat]["g1"]
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        th=X[4]
        Me=X[5]
        A=np.pi*d*d/4
        g=Me*0.3048*10**self.C1*(61024*b*A)**self.a1*(15432.4*m)**self.b1*(1/np.cos(th))**self.g1-v
        #g=0.3048*10**C1*(61024*b*A)**a1*(2.2046*m)**b1*(1/np.cos(th))**g1-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        th=X[4]
        Me=X[5]
        C1=self.C1
        b1=self.b1
        a1=self.a1
        g1=self.g1
        dGdX[0]=Me*0.3048*10**C1*a1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1/b
        dGdX[1] =Me*0.6096*10**C1*a1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1/d
        dGdX[2] =Me*0.3048*10**C1*b1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1/m
        dGdX[3] =-1
        dGdX[4] =Me*0.3048*10**C1*g1*(15432.4*m)**b1*(47928.1375231659*b*d**2)**a1*(1/np.cos(th))**g1*np.sin(th)/np.cos(th)
        dGdX[5]=0.3048*10**C1*(61024*b*np.pi*d*d/4)**a1*(15432.4*m)**b1*(1/np.cos(th))**g1
        return dGdX
################################
#             BRL            #
################################
class BRL(penMed):
    """
    ---  Ballistic Research Laboratories (BRL) model (1968)  ---
    GG Corbett and SR Reid.,Quasi-static and dynamic local
    loading of monolithic simply supported steel plate.
    International journal of impact engineering, Vol. 13, No. 3, pp.
    423–441, 1993.
        ---
    ***variables***
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    Me  Model error
    ***constants***
    Limp    length of impactor
    Lsh     unsupported shield panel span
    Su      ultimate tensile strength of shield material
    """
    def __init__(self):
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','v','Me']
        self.const=['Limp','Lsh','Su']
        title='BRL Formula'
        val_range={
            'v_bl':[57,270],
            'Limp/d':[1.25,8],
            'b/d':[0.1,1.0],
            'Lsh/d':[8,35],
            'Su':[315,500]
        }
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
        super().SaveConst(self.const)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Limp=data['Limp']['mean']
        Lsh=data['Lsh']['mean']
        Su=data['Su']['mean']
        a7=5.37
        v_bl=a7*1e4*(b*d)**0.75/m**0.5
        if self.i_Valid:
            super().check('v_bl',v_bl)
            super().check('Limp/d',Limp/d)
            super().check('b/d',b/d)
            super().check('Lsh/d',Lsh/d)
            super().check('Su',Su)
        else:
            ii=0
            if super().check_c('v_bl',v_bl)!=True: ii+=1
            if super().check_c('Limp/d',Limp/d)!=True: ii+=1
            if super().check_c('b/d',b/d)!=True: ii+=1
            if super().check_c('Lsh/d',Lsh/d)!=True: ii+=1
            if super().check_c('Su',Su)!=True: ii+=1
            return ii

    def gcalc(self,X):
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Me=X[4]
        a7=5.37
        g=Me*a7*1e4*(b*d)**0.75/m**0.5-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Me=X[4]
        a7=5.37
        dGdX[0]= Me*7500.0*a7*(b*d)**0.75/(b*m**0.5)
        dGdX[1]= Me*7500.0*a7*(b*d)**0.75/(d*m**0.5)
        dGdX[2]= Me*(-5000.0)*a7*(b*d)**0.75/m**1.5
        dGdX[3] =-1
        dGdX[4] = 10000.0*a7*(b*d)**0.75/m**0.5
        return dGdX
################################
#             DeMarre        #
################################
class DeMarre(penMed):
    """
    ---  De Marre formula (Herrmann and Jones,1961)  ---
    Walter Herrmann and Arfon H Jones.
    Survey of hypervelocity impact information. ASRL, 1961.
    ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    """
    def __init__(self):
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','v','Me']
        self.const=[]
        title='De Marre Formula'
        val_range={
            'v_bl':[200,900],
            'm':[0.1,50]
        }
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        v_bl=0.4311e5*d**0.75*b**0.7/m**0.5
        if self.i_Valid:
            super().check('v_bl',v_bl)
            super().check('m',m)
        else:
            ii=0
            if super().check_c('v_bl',v_bl)!=True: ii+=1
            if super().check_c('m',m) !=True: ii+=1
            return ii           
    def gcalc(self,X):
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Me=X[4]
        g=Me*0.431e5*d**0.75*b**0.7/m**0.5-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Me=X[4]
        dGdX[0]=Me*30170.0*d**0.75/(b**0.3*m**0.5)
        dGdX[1] =Me*32325.0*b**0.7/(d**0.25*m**0.5)
        dGdX[2] =Me*(-21550.0)*b**0.7*d**0.75/m**1.5
        dGdX[3] =-1
        dGdX[4] =0.431e5*d**0.75*b**0.7/m**0.5
        return dGdX
################################
#             Jowett         #
################################
class Jowett(penMed):
    """
    ---  Jowett Formula (1986)  ---
    GG Corbett, SR Reid, and W Johnson. 
    Impact loading of plates and shells by free-flying
    projectiles: a review. International Journal 
    of Impact Engineering, Vol. 18, No. 2, pp.141–230, 1996.
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
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','Su','v','Me']
        self.const=['Lsh','Limp']
        title='Jowett Formula'
        val_range={
            'vbl':[40,200],
            'Su':[315,483],
            'Limp/d':[2,8],
            'b/d':[0.1,0.64]
        }
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
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
        if self.i_Valid:
            super().check('vbl',vbl)
            super().check('Su',Su)
            super().check('Limp/d',Limp/d)
            super().check('b/d',b/d)
        else:
            ii=0
            if super().check_c('vbl',vbl)!=True: ii+=1
            if super().check_c('Su',Su)!=True: ii+=1
            if super().check_c('Limp/d',Limp/d)!=True: ii+=1
            if super().check_c('b/d',b/d)!=True: ii+=1
            return ii
        
    def gcalc(self,X):
        df=super().gDict()
        b=X[0]
        d=X[1]
        m=X[2]
        Su=X[3]
        v=X[4]
        Me=X[5]
        ratio=b/d
        self.Lsh=df['Lsh']['mean']
        if self.Lsh/d<=12:
            omg=(self.Lsh/d)**0.305
        else:
            #omg=12.0
            omg=2.14
        self.omg=omg
        ratio=b/d
        vbl=0
        if ratio>0.1 and ratio <0.25:
            vbl=Me*1.62*omg*d*np.sqrt(Su*d/m)*(b/d)**0.87
        if ratio>=0.25 and ratio<0.64:
            vbl=Me*0.87*omg*d*np.sqrt(Su*d/m)*(b/d)**0.42
        g=vbl-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        Su=X[3]
        v=X[4]
        Me=X[5]
        ratio=b/d
        omg=self.omg
        if ratio>0.1 and ratio <0.25:
            dGdX[0]= Me*1.4094*d*omg*(b/d)**0.87*np.sqrt(Su*d/m)/b
            dGdX[1]= Me*1.0206*omg*(b/d)**0.87*np.sqrt(Su*d/m)
            dGdX[2]= Me*(-0.81)*d*omg*(b/d)**0.87*np.sqrt(Su*d/m)/m
            dGdX[3]= Me*0.81*d*omg*(b/d)**0.87*np.sqrt(Su*d/m)/Su
            dGdX[4]= -1
            dGdX[5]= 1.62*omg*d*np.sqrt(Su*d/m)*(b/d)**0.87
        if ratio>=0.25 and ratio<0.64:
            dGdX[0]= Me*0.3654*d*omg*(b/d)**0.42*np.sqrt(Su*d/m)/b
            dGdX[1]= Me*0.9396*omg*(b/d)**0.42*np.sqrt(Su*d/m)
            dGdX[2]= Me*(-0.435)*d*omg*(b/d)**0.42*np.sqrt(Su*d/m)/m
            dGdX[3]= Me*0.435*d*omg*(b/d)**0.42*np.sqrt(Su*d/m)/Su
            dGdX[4]= -1
            dGdX[5]= Me*0.87*omg*d*np.sqrt(Su*d/m)*(b/d)**0.42
        return dGdX
################################
#             Lambert        #
################################
class Lambert(penMed):
    """
    ---  Lambert and Jonas Approximation (1976)  ---
    JP Lambert and GH Jonas.
    Towards standardization in terminal ballistics testing:
    velocity representation.
    USA Ballistic Research Laboratories, Aberdeen Proving Ground, MD,
    Technical Report, No. 1852, 1976.
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
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','th','v','Limp','Me']
        self.const=['ro_imp','Material']
        title='Lambert and Jonas Approximation'
        val_range={
            'm':[0.0005,3.63],
            'd':[0.002,0.05],
            'Limp/d':[4,30],
            'b':[0.006,0.15],
            'th':[0,60/360],
            'ro_imp':[7800,19000]
        }
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        th=data['th']['mean']
        Limp=data['Limp']['mean']
        ro_imp=data['ro_imp']['mean']
        mat=data['Material']
        #self.setMaterial(mat)
        if self.i_Valid:
            super().check('b',b)
            super().check('d',d)
            super().check('m',m)
            super().check('th',th)
            super().check('Limp/d',Limp/d)
            super().check('ro_imp',ro_imp)
        else:
            ii=0
            if super().check_c('b',b)!=True: ii+=1
            if super().check_c('d',d)!=True: ii+=1
            if super().check_c('m',m)!=True: ii+=1
            if super().check_c('th',th)!=True: ii+=1
            if super().check_c('Limp/d',Limp/d)!=True: ii+=1
            if super().check_c('ro_imp',ro_imp)!=True: ii+=1
            return ii
    def MatList(self):
        return ['aluminum','RHA']
    def gcalc(self,X):
        df=super().gDict()
        mat=df['Material']
        if mat=='aluminum':
            a10=1750
        if mat=='RHA':
            a10=4000
        self.a10=a10
        b=X[0]
        d=X[1]
        m=X[2]
        th=X[3]
        v=X[4]
        Limp=X[5]
        Me=X[6]
        z=(b/d)*(1/np.cos(th))**0.75
        f=z+np.exp(z)-1
        g=Me*31.62*a10*(Limp/d)**0.15*np.sqrt(f*d**3/m)-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        th=X[3]
        v=X[4]
        Limp=X[5]
        Me=X[6]
        a10=self.a10
        z=(b/d)*(1/np.cos(th))**0.75
        f=z+np.exp(z)-1
        dGdX[0]= Me*15.81*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)*((1/np.cos(th))**0.75*np.exp(b*(1/np.cos(th))**0.75/d)/d + (1/np.cos(th))**0.75/d)/(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)
        dGdX[1]= Me*(-4.743)*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)/d + 31.62*a10*m*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)*(d**3*(-b*(1/np.cos(th))**0.75*np.exp(b*(1/np.cos(th))**0.75/d)/d**2 - b*(1/np.cos(th))**0.75/d**2)/(2*m) + 3*d**2*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/(2*m))/(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1))
        dGdX[2]= Me*(-15.81)*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)/m
        dGdX[3]= Me*15.81*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)*(0.75*b*(1/np.cos(th))**0.75*np.exp(b*(1/np.cos(th))**0.75/d)*np.sin(th)/(d*np.cos(th)) + 0.75*b*(1/np.cos(th))**0.75*np.sin(th)/(d*np.cos(th)))/(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)
        dGdX[4]= -1
        dGdX[5]= Me*4.743*a10*(Limp/d)**0.15*np.sqrt(d**3*(b*(1/np.cos(th))**0.75/d + np.exp(b*(1/np.cos(th))**0.75/d) - 1)/m)/Limp
        dGdX[6]=31.62*a10*(Limp/d)**0.15*np.sqrt(f*d**3/m)
        return dGdX
################################
#             Neilson        #
################################
class Neilson(penMed):
    """
    ---  Neilson Formula (1985)  ---
    AJ Neilson. Empirical equations for the perforation
    of mild steel plates. 
    International Journal of Impact Engineering, Vol. 3, No. 2, pp. 137–142, 1985.
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
    shape 'flat' or 'hemispherical' (ノーズ形状)
    """
    def __init__(self):
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','Su','Lsh','v','Me']
        self.const=['Limp','shape']
        title='Neilson Formula'
        val_range={
            'b/d':[0.14,0.64],
            'Lsh/d':[4,22],
            'Limp/d':[13,1e4],
        }
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        Limp=data['Limp']['mean']
        shape=data['shape']
        #self.SetNose(Lsh,d,shape)
        if self.i_Valid:
            super().check('b/d',b/d)
            #super().check('Lsh/d',Lsh/d)
            super().check('Limp/d',Limp/d)
        else:
            ii=0
            if super().check_c('b/d',b/d)!=True:ii+=1
            if super().check_c('Limp/d',Limp/d)!=True:ii+=1
            return ii
    def gcalc(self,X):
        df=super().gDict()
        Lsh=df['Lsh']['mean']
        d=df['d']['mean']
        shape=df['shape']
        if shape=='flat':
            if Lsh/d > 4.0 and Lsh/d<22.0:
                a12=1.67
            if Lsh/d >=22.0:
                a12=4.26
        a12=1.67 #上の条件を無視
        if shape=='hemispherical':
            a12=4.24
        self.a12=a12          
        b=X[0]
        d=X[1]
        m=X[2]
        Su=X[3]
        Lsh=X[4]
        v=X[5]
        Me=X[6]
        g=Me*a12*d*np.sqrt(Su*d/m)*(b/d)**0.85*(Lsh/d)**0.3-v
        return g
    def dGdXcalc(self,X):
        a12=self.a12
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        Su=X[3]
        Lsh=X[4]
        v=X[5]
        Me=X[6]
        dGdX[0]= Me*0.85*a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/b
        dGdX[1]= Me*0.35*a12*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)
        dGdX[2]= Me*(-a12)*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/(2*m)
        dGdX[3]= Me*a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/(2*Su)
        dGdX[4]= Me*0.3*a12*d*(Lsh/d)**0.3*(b/d)**0.85*np.sqrt(Su*d/m)/Lsh
        dGdX[5]= -1
        dGdX[6]=a12*d*np.sqrt(Su*d/m)*(b/d)**0.85*(Lsh/d)**0.3
        return dGdX
################################
#             Ohte           #
################################
class Ohte(penMed):
    """
    ---  Ohte et al. Formula (Ohte et al., 1982)  ---
    Satoshi OHTE, Hiroyasu YOSHIZAWA, Norimasa CHIBA, and Shigeru SHIDA. 
    Impact strength of steel plates struck by projectiles:
    evaluation formula for critical fracture energy of steel plate.
    Bulletin of JSME, Vol. 25, No. 206, pp. 1226–1231, 1982.
        ---
    b   thickness of a shield
    d   maximum diameter of impactor
    m   initial mass of the impactor
    v   velocity of impactor
    """
    def __init__(self):
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','v','Me']
        self.const=['Lsh','Su']
        title='Ohte et al. Formula'
        val_range={
            'v_bl':[25,180],
            'm':[3,50],
            'Su':[490,637],
            'b':[7,38],
            'Lsh/b':[39,1e4],
            'd':[39,1e4]
        }
        super().SaveTitle(title)
        super().SaveRange(val_range)
        super().SaveVariable(self.variable)
    def Validation(self,data):
        b=data['b']['mean']
        d=data['d']['mean']
        m=data['m']['mean']
        Su=data['Su']['mean']
        Lsh=data['Lsh']['mean']
        v_bl=7.67e4*(b*d)**0.75/m**0.5
        if self.i_Valid:
            super().check('v_bl',v_bl)
            super().check('m',m)
            super().check('Su',Su)
            super().check('b',b)
            super().check('Lsh/b',Lsh/b)
            super().check('d',d)
        else:
            ii=0
            if super().check_c('v_bl',v_bl)!=True:ii+=1
            if super().check_c('m',m)!=True:ii+=1
            if super().check_c('Su',Su)!=True:ii+=1
            if super().check_c('b',b)!=True:ii+=1
            if super().check_c('Lsh/b',Lsh/b)!=True:ii+=1
            if super().check_c('d',d)!=True:ii+=1
            return ii

    def gcalc(self,X):
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Me=X[4]
        g=Me*7.67e4*(b*d)**0.75/m**0.5-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Me=X[4]
        dGdX[0]= Me*57525.0*(b*d)**0.75/(b*m**0.5)
        dGdX[1]= Me*57525.0*(b*d)**0.75/(d*m**0.5)
        dGdX[2]= Me*(-38350.0)*(b*d)**0.75/m**1.5
        dGdX[3]= -1
        dGdX[4]=7.67e4*(b*d)**0.75/m**0.5
        return dGdX
################################
#             SRI            #
################################
class SRI(penMed):
    """
    ---  Stanford Research Institute (SRI) correlation (1963)  ---
    GG Corbett, SR Reid, and W Johnson. 
    Impact loading of plates and shells by free-flying
    projectiles: a review. 
    International Journal of Impact Engineering, Vol. 18, No. 2, pp.
    141–230, 1996.
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
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','v','Lsh','Su','Me']
        self.const=['Limp']
        title='SRI Formula'
        val_range={
            'v_bl':[21,122],
            'b/d':[0.1,0.6],
            'Lsh/d':[5,8],
            'b/Lsh':[0.002,0.05],
            'Lsh/b':[0.0,100],
            'Limp/d':[5,8]
        }
        super().SaveTitle(title)
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
        if self.i_Valid:
            super().check('v_bl',v_bl)
            super().check('b/d',b/d)
            super().check('Lsh/d',Lsh/d)
            super().check('b/Lsh',b/Lsh)
            super().check('Lsh/b',Lsh/b)
            super().check('Limp/d',Limp/d)
        else:
            ii=0
            if super().check_c('v_bl',v_bl)!=True:ii+=1
            if super().check_c('b/d',b/d)!=True:ii+=1
            if super().check_c('Lsh/d',Lsh/d)!=True:ii+=1
            if super().check_c('b/Lsh',b/Lsh)!=True:ii+=1
            if super().check_c('Lsh/b',Lsh/b)!=True:ii+=1
            if super().check_c('Limp/d',Limp/d) !=True:ii+=1
            return ii           
    def gcalc(self,X):
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Su=X[4]
        Lsh=X[5]
        Me=X[6]
        a6=0.4
        g=Me*a6*b*np.sqrt(Su*d/m*(42.7+Lsh/b))-v
        return g
    def dGdXcalc(self,X):
        a6=0.4
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        v=X[3]
        Su=X[4]
        Lsh=X[5]
        Me=X[6]
        dGdX[0]= -Me*Lsh*a6*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*b*(Lsh/b + 42.7)) + a6*np.sqrt(Su*d*(Lsh/b + 42.7)/m)
        dGdX[1]= Me*a6*b*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*d)
        dGdX[2]= -Me*a6*b*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*m)
        dGdX[3]= -1
        dGdX[4]= Me*a6*b*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*Su)
        dGdX[5]= Me*a6*np.sqrt(Su*d*(Lsh/b + 42.7)/m)/(2*(Lsh/b + 42.7))
        dGdX[6]=a6*b*np.sqrt(Su*d/m*(42.7+Lsh/b))
        return dGdX
################################
#             SwRI           #
################################
class SwRI(penMed):
    """
    ---  Southwest Research Institute (SwRI) model ---
    W.E. Baker, J.J. Kulesz, P.S. Westine, P.A. Cox, and J.S. Wilbeck.
    A manual for the prediction of blast and fragment loading on structures.
    Report DOE/TIC-11268. United States Department of Energy,
    Albuquerque Operation Office, Amarillo Area Office, Amarillo,TX, 1980.
    ---
    b   thickness of a shield
    m   initial mass of the impactor
    v   velocity of impactor
    th  angle between a normal vector to a shield surface and the direction of impactor
    fragment  :'Standard' or 'Alternative'
    """
    def __init__(self):
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','m','v','th','Me']
        title='Southwest Research Institute (SwRI) model'
        self.const=['fragment']
        super().SaveTitle(title)
        super().SaveRange('Validation process is not defined.')
        super().SaveVariable(self.variable)
    def gTable(self,df):#このclass用の特別処理
        tab={"0":{"b1":1414,"b2":0.295,"b3":0.910},
        "1":{"b1":1936,"b2":0.096,"b3":1.310},
        "2":{"b1":2039,"b2":0.064,"b3":0.430}}
        m=df['m']['mean']
        b=df['b']['mean']
        if df['fragment'] == 'Standard':
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
        return S,b1,b2,b3
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
        if self.i_Valid:
            print('Validation process is not defined.')
        else:
            return 0
    def gcalc(self,X):
        df=super().gDict()
        S,b1,b2,b3=self.gTable(df)
        b=X[0]
        m=X[1]
        v=X[2]
        th=X[3]
        Me=X[4]
        g=Me*0.205*b1/np.sqrt(m)*S**b2*(39.37*b/np.cos(th))**b3-v
        return g
    def dGdXcalc(self,X):
        df=super().gDict()
        S,b1,b2,b3=self.gTable(df)
        dGdX=[0]*len(X)
        b=X[0]
        m=X[1]
        v=X[2]
        th=X[3]
        Me=X[4]
        dGdX[0]= Me*0.205*S**b2*b1*b3*(39.37*b/np.cos(th))**b3/(b*np.sqrt(m))
        dGdX[1]= Me*(-0.1025)*S**b2*b1*(39.37*b/np.cos(th))**b3/m**(3/2)
        dGdX[2]= -1
        dGdX[3]= Me*0.205*S**b2*b1*b3*(39.37*b/np.cos(th))**b3*np.sin(th)/(np.sqrt(m)*np.cos(th))
        dGdX[4]=0.205*b1/np.sqrt(m)*S**b2*(39.37*b/np.cos(th))**b3
        return dGdX
################################
#             WenJones       #
################################
class WenJones(penMed):
    """
    ---   Wen and Jones Formula (1992)   ---
    T. Børvik, M. Langseth, O.S. Hopperstad, and K.A. Malo.
    Empirical equations for ballistic penetration of metal plates. 
    Fortifikatorisk Notat No.260/98. The Norwegian Defence
    Construction Service, Central Staff - Technical Division, Oslo, Norway, 1998.
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
        self.i_Valid=True  #Validation結果を出力するときTrue
        self.variable=['b','d','m','Sy','Lsh','v','Me']
        self.const=['Su']
        title='Wen Jones Formula'
        val_range={
            'vbl':[0,20],
            'Su':[340,440],
            'Lsh/d':[40,40],
            'Lsh/b':[25,100],
            'b/d':[0.4,1.6]
        }
        super().SaveTitle(title)
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
        if self.i_Valid:
            super().check('vbl',vbl)
            super().check('Su',Su)
            super().check('Lsh/d',Lsh/d)
            super().check('Lsh/b',Lsh/b)
            super().check('b/d',b/d)
        else:
            ii=0
            if super().check_c('vbl',vbl)!=True:ii+=1
            if super().check_c('Su',Su)!=True:ii+=1
            if super().check_c('Lsh/d',Lsh/d)!=True:ii+=1
            if super().check_c('Lsh/b',Lsh/b)!=True:ii+=1
            if super().check_c('b/d',b/d) !=True:ii+=1
            return ii          

    def gcalc(self,X):
        b=X[0]
        d=X[1]
        m=X[2]
        Sy=X[3]
        Lsh=X[4]
        v=X[5]
        Me=X[6]
        g=Me*2*d*np.sqrt(Sy*d/m*(0.25*np.pi*(b/d)**2+(b/d)**1.47*(Lsh/d)**0.21))-v
        return g
    def dGdXcalc(self,X):
        dGdX=[0]*len(X)
        b=X[0]
        d=X[1]
        m=X[2]
        Sy=X[3]
        Lsh=X[4]
        v=X[5]
        Me=X[6]
        dGdX[0]= Me*d*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)*(1.5707963267949*b/d**2 + 1.47*(Lsh/d)**0.21*(b/d)**1.47/b)/(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)
        dGdX[1]= Me*2*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m) + 2*m*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)*(Sy*d*(-1.5707963267949*b**2/d**3 - 1.68*(Lsh/d)**0.21*(b/d)**1.47/d)/(2*m) + Sy*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/(2*m))/(Sy*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47))
        dGdX[2]= Me*(-d)*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)/m
        dGdX[3]= Me*d*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)/Sy
        dGdX[4]= Me*0.21*d*(Lsh/d)**0.21*(b/d)**1.47*np.sqrt(Sy*d*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47)/m)/(Lsh*(0.785398163397448*b**2/d**2 + (Lsh/d)**0.21*(b/d)**1.47))
        dGdX[5]= -1
        dGdX[6]=2*d*np.sqrt(Sy*d/m*(0.25*np.pi*(b/d)**2+(b/d)**1.47*(Lsh/d)**0.21))
        return dGdX