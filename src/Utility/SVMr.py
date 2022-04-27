from sklearn.svm import SVC
import numpy as np
import math
class SVMr():
    """
    *** SVMを用いた信頼性解析管理クラス ***
    関連文献:
    酒井信介、木原重光、「サポートベクタマシンを用いたプラント保全データの信頼性評価」
                   、圧力技術、Vol.59、No.3，pp.129-139(2021).
        例題プログラム:データtrain_xとそのクラスtrain_t1が与えられているとき        
        svmr=SVMr(train_x,train_t1) #インスタンスの生成
        svmr.SVM() #SVM解析
        dp,beta=svmr.SV_RF() #提案アルゴリズムによる設計点探索、設計点、信頼性指標の出力
        alpha,deriv=svmr.sv_alpha(dp)#設計点での感度
    """
    def __init__(self,X,y):
        self.X=X
        self.y=y
        self.gamma=0.1
        self.kernel='rbf'
        self.C=10
    def SetGamma(self,x):
        self.gamma=x
    def SetKernel(self,kernel):
        self.kernel=kernel
    def SetC(self,x):
        self.C=x
    def SVM(self):
        self.svm=SVC(kernel=self.kernel,C=self.C,gamma=self.gamma,probability=True).fit(self.X,self.y)
        self.sv=svm.support_vectors_
    def SV_RF(self):
        """
        目的:全SV点を出発点とするRFを実施し、最小βを与える点を設計点と確定
        入力:
            sv  SVM解析の結果得られるsv点
            svm     SVM解析の出力
            gamma   SVM解析の際に用いたγ値
        出力:dp,beta
            Dp      原点から最短の設計点
            Beta　  β値
        """
        beta=np.zeros(len(self.sv))
        for i in range(len(self.sv)):
            dp,beta[i]=self.RackwitzFiessler(self.sv[i])
        start=self.sv[np.argmin(beta)]
        Dp,Beta=self.RackwitzFiessler(start)
        return Dp,Beta
    def RackwitzFiessler(self,start):
        """
        目的:初期点startからスタートし，Rackwitz Fiessler法により設計点を求める
        入力:
            start       初期点
            svm     SVM解析の出力
            gamma   SVM解析の際に用いたγ値
            nmax    繰り返し数上限
            eps     収束規準
            b0      beta値初期値
        出力:dp,beta
            dp      設計点
            beta    信頼性指標
        """
        x=start
        nmax=100; eps=0.001; b0=10.0
        for i in range(nmax):
            alpha,nabla=self.sv_alpha(x)
            if(np.isnan(alpha[0])):
                beta=1e4
                break
            x2=1/np.dot(nabla,nabla)*(np.dot(nabla,x)-self.g(x))*nabla.T
            x=x2
            beta=np.linalg.norm(x)
            de=abs((beta-b0)/b0)
            b0=beta
            if de<eps:
                break
        dp=x
        return dp,beta
    def sv_alpha(self,x0):
        """
        目的:座標点x0の感度係数と，微分値を計算する
        入力:
            svm SVM解析の出力
            gamma   SVM解析の際に用いたγ値
            x0      評価点の座標
        出力:alpha,deriv
            alpha   感度ベクトル
            deriv   微分値ベクトル
        """
        #derivative
        var_num=len(x0)
        coef=self.svm.dual_coef_
        sv=self.svm.support_vectors_
        b=self.svm.intercept_
        # value of kernel function
        nsv=len(sv)
        kernel=np.zeros(nsv)
        #x0=np.array(dp)
        for i in range(nsv):
            ee=-self.gamma*np.linalg.norm(x0-sv[i])**2
            kernel[i]=math.exp(ee)
        #derivative
        deriv=np.zeros(var_num)
        for j in range(var_num):
            ar_d=(-2.*self.gamma)*coef*(x0[j]-sv[:,j])
            deriv[j]=np.dot(ar_d,kernel.T)
        alpha=deriv/np.linalg.norm(deriv)
        # alpha: 感度係数
        # deriv: 各変数方向の微分値
        return alpha,deriv
    def g(self,x):    # svmでのgのsurrogate関数
        """
        目的:SVMの解析結果から，x点のdecision_functionの値を返す
        入力:
            x   評価点
            svm SVM解析の出力
        出力:
            decisio_functionの値
        """
        var_num=len(x)
        return self.svm.decision_function(x.reshape(1,var_num))[0]