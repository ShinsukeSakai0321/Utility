import tkinter as tk
from tkinter import ttk
#from Utility import Penetration as pen
import Penetration as pen
import pandas as pd
import pickle
class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title('Penetration analysis system')
        self.method=[
            'AlyLi',
            'BRL',
            'DeMarre',
            'Jowett',
            'Lambert',
            'Neilson',
            'Ohte',
            'SRI',
            'SwRI',
            'THOR',
            'WenJones'  
        ]
        lbl = tk.Label(text='Formula')
        lbl.place(x=20,y=10)
        # チェックボックスON/OFFの状態

        self.var_item = tk.IntVar()
        self.var_prob=tk.IntVar()
        #手法の設定
        for i in range(len(self.method)):
            btn=tk.Radiobutton(self.master,value=i, variable=self.var_item,text=self.method[i],command=self.change_selected_item) 
            btn.place(x=10, y=30 + (i * 24))
        btn_v=tk.Radiobutton(self.master,value=0,variable=self.var_prob,text='Evaluate v_lim',command=self.change_selected_item)
        btn_p=tk.Radiobutton(self.master,value=1,variable=self.var_prob,text='Probabilistic analysis',command=self.change_selected_item)
        #self.var=['b', 'd', 'm', 'Su', 'Sy','Lsh', 'th', 'Limp', 'v','ro_imp','a10','shape','fragment','Material']
        self.var=['b', 'd', 'm', 'Su', 'Sy','Lsh', 'th', 'Limp', 'v','ro_imp','shape','fragment','Material']
        lbl_b = tk.Label(text='b'); self.lbl_b=tk.Label(text='*')
        lbl_d = tk.Label(text='d'); self.lbl_d=tk.Label(text='*')
        lbl_m = tk.Label(text='m'); self.lbl_m=tk.Label(text='*')
        lbl_Su = tk.Label(text='Su'); self.lbl_Su=tk.Label(text='*')
        lbl_Sy = tk.Label(text='Sy'); self.lbl_Sy=tk.Label(text='*')
        lbl_Lsh = tk.Label(text='Lsh'); self.lbl_Lsh=tk.Label(text='*')
        lbl_th = tk.Label(text='th'); self.lbl_th=tk.Label(text='*')
        lbl_Limp = tk.Label(text='Limp'); self.lbl_Limp=tk.Label(text='*')
        lbl_v = tk.Label(text='v'); self.lbl_v=tk.Label(text='*')
        lbl_ro_imp = tk.Label(text='ro_imp'); self.lbl_ro_imp=tk.Label(text='*')
        lbl_a10 = tk.Label(text='a10'); self.lbl_a10=tk.Label(text='*')
        lbl_shape=tk.Label(text='Shape');self.lbl_shape=tk.Label(text='*')
        lbl_frag=tk.Label(text='Fragment');self.lbl_frag=tk.Label(text='*')
        lbl_Material=tk.Label(text='Material');self.lbl_Material=tk.Label(text='*')
        lbl_const=tk.Label(text='Const.')
        lbl_mean=tk.Label(text='Mean')
        lbl_cov=tk.Label(text='COV')
        en=15
        self.txt_b_m = tk.Entry(width=en)
        self.txt_d_m = tk.Entry(width=en)
        self.txt_m_m = tk.Entry(width=en)
        self.txt_Su_m = tk.Entry(width=en)
        self.txt_Sy_m = tk.Entry(width=en)
        self.txt_Lsh_m = tk.Entry(width=en)
        self.txt_th_m = tk.Entry(width=en)
        self.txt_Limp_m = tk.Entry(width=en)
        self.txt_v_m = tk.Entry(width=en)
        self.txt_b_c = tk.Entry(width=en)
        self.txt_d_c = tk.Entry(width=en)
        self.txt_m_c = tk.Entry(width=en)
        self.txt_Su_c = tk.Entry(width=en)
        self.txt_Sy_c = tk.Entry(width=en)
        self.txt_Lsh_c = tk.Entry(width=en)
        self.txt_th_c = tk.Entry(width=en)
        self.txt_Limp_c = tk.Entry(width=en)
        self.txt_v_c = tk.Entry(width=en)
        self.txt_ro_imp=tk.Entry(width=en)
        
        #self.txt_a10=tk.Entry(width=en*2)
        #self.txt_shape=tk.Entry(width=en*2)
        self.combo_shape=ttk.Combobox(width=en,values=['none'])
        #self.txt_frag=tk.Entry(width=en*2)
        self.combo_frag=ttk.Combobox(width=en,values=['none'])
        self.combo_mat=ttk.Combobox(width=en,values=['none'])
        self.txt_Text=tk.Text(height=25,width=70)
        
        lbl_title=tk.Label(text='Title')
        self.txt_title=tk.Entry(width=en*5)
        button_calc = tk.Button(self.master, text = "Calc", command = self.Calc_click)
        button_load= tk.Button(self.master, text = "Load", command = self.Load_click)
        self.txt_load=tk.Entry(width=en)
        button_save= tk.Button(self.master, text = "Save", command = self.Save_click)
        button_exit=tk.Button(self.master,text="Exit",command=self.master.destroy)
        self.txt_save=tk.Entry(width=en)
        i=0; xx=150; x_ast=xx-10; yy=30; dy_ast=3; dx0=50; dx=100
        button_exit.place(x=30,y=30 + (len(self.method) * 24)+70)
        self.txt_Text.place(x=xx+260,y=yy)
        lbl_const.place(x=xx-30,y=yy-24)
        lbl_mean.place(x=xx+dx0,y=yy-24)
        lbl_cov.place(x=xx+dx0+dx,y=yy-24)
        lbl_b.place(x=xx,y=yy);  self.lbl_b.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_b_m.place(x=xx+dx0,y=yy)
        self.txt_b_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_d.place(x=xx,y=yy); self.lbl_d.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_d_m.place(x=xx+dx0,y=yy)
        self.txt_d_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_m.place(x=xx,y=yy); self.lbl_m.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_m_m.place(x=xx+dx0,y=yy)
        self.txt_m_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_Su.place(x=xx,y=yy); self.lbl_Su.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_Su_m.place(x=xx+dx0,y=yy)
        self.txt_Su_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_Sy.place(x=xx,y=yy); self.lbl_Sy.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_Sy_m.place(x=xx+dx0,y=yy)
        self.txt_Sy_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_Lsh.place(x=xx,y=yy); self.lbl_Lsh.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_Lsh_m.place(x=xx+dx0,y=yy)
        self.txt_Lsh_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_th.place(x=xx,y=yy); self.lbl_th.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_th_m.place(x=xx+dx0,y=yy)
        self.txt_th_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_Limp.place(x=xx,y=yy); self.lbl_Limp.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_Limp_m.place(x=xx+dx0,y=yy)
        self.txt_Limp_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_v.place(x=xx,y=yy); self.lbl_v.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_v_m.place(x=xx+dx0,y=yy)
        self.txt_v_c.place(x=xx+dx0+dx,y=yy)
        yy+=24;lbl_ro_imp.place(x=xx,y=yy); self.lbl_ro_imp.place(x=x_ast,y=yy+dy_ast);i+=1
        self.txt_ro_imp.place(x=xx+dx0,y=yy)
        #yy+=24;lbl_a10.place(x=xx,y=yy); self.lbl_a10.place(x=x_ast,y=yy+dy_ast);i+=1
        #self.txt_a10.place(x=xx+dx0,y=yy)
        yy+=24;lbl_shape.place(x=xx,y=yy); self.lbl_shape.place(x=x_ast,y=yy+dy_ast);i+=1
        self.combo_shape.place(x=xx+dx0,y=yy)
        yy+=24;lbl_frag.place(x=xx,y=yy); self.lbl_frag.place(x=x_ast,y=yy+dy_ast);i+=1
        self.combo_frag.place(x=xx+dx0,y=yy)
        yy+=24;lbl_Material.place(x=xx,y=yy); self.lbl_Material.place(x=x_ast,y=yy+dy_ast);i+=1
        self.combo_mat.place(x=xx+dx0,y=yy)
        yy+=48;btn_v.place(x=xx,y=yy)
        btn_p.place(x=xx+dx,y=yy)
        yy+=48; lbl_title.place(x=xx,y=yy); self.txt_title.place(x=xx+dx0,y=yy)
        yy+=48; button_calc.place(x=xx,y=yy)
        yy+=48; button_load.place(x=xx,y=yy); self.txt_load.place(x=xx+dx0,y=yy)
        yy+=48; button_save.place(x=xx,y=yy); self.txt_save.place(x=xx+dx0,y=yy)

        self.change_selected_item()
    def Save_click(self):
        df=self.MakeDict()
        fname=self.txt_save.get()
        with open(fname, mode='wb') as f:
            pickle.dump(df,f)
    def Load_click(self):
        fname=self.txt_load.get()
        with open(fname,'rb') as f:
            df = pickle.load(f)
        self.toDict(df)
    def Calc_click(self):
        df=self.MakeDict()
        if 'Material' in df.keys():
            self.formula.setMaterial(df['Material'])
        print('--------------',self.formula.title,'--------------')
        self.formula.Validation(df)
        if self.iProb==1:
            print('*** Probabilistic analysis ***')
            print('variable=',self.formula.variable)
            print('value=',self.formula.Gcheck(df))#確率変数の平均値に対するg値評価
            self.formula.Reliability(df)#信頼性評価
            print('beta=',self.formula.GetBeta())#信頼性指標の出力
            print('Alpha=',self.formula.GetAlpha())#感度の出力
            print('Pf=',self.formula.GetPOF())#破損確率の出力
        else:
            print('*** Analysis of Balistic Limit Velocity ***')
            df['v']['mean']=0
            print('Vbl=',self.formula.Gcheck(df))
        return
    def MakeDict(self):
        df={}
        df=self.dfAdd(df)
        df=self.dfAddC(df)
        df['Title']=self.txt_title.get()
        return df
    def Clear(self):
        self.txt_title.delete(0,tk.END)
        self.txt_b_m.delete(0,tk.END)
        self.txt_b_c.delete(0,tk.END)
        self.txt_d_m.delete(0,tk.END)
        self.txt_d_c.delete(0,tk.END)
        self.txt_m_m.delete(0,tk.END)
        self.txt_m_c.delete(0,tk.END)
        self.txt_Su_m.delete(0,tk.END)
        self.txt_Su_c.delete(0,tk.END)
        self.txt_Sy_m.delete(0,tk.END)
        self.txt_Sy_c.delete(0,tk.END)
        self.txt_Lsh_m.delete(0,tk.END)
        self.txt_Lsh_c.delete(0,tk.END)
        self.txt_th_m.delete(0,tk.END)
        self.txt_th_c.delete(0,tk.END)
        self.txt_Limp_m.delete(0,tk.END)
        self.txt_Limp_c.delete(0,tk.END)
        self.txt_v_m.delete(0,tk.END)
        self.txt_v_c.delete(0,tk.END)
        self.txt_ro_imp.delete(0,tk.END)
        #self.txt_a10.delete(0,tk.END)
        #self.combo_shape.delete(0,tk.END)
        #self.combo_frag.delete(0,tk.END)

    def toDict(self,df):
        self.Clear()
        self.txt_title.insert(0,df['Title'])
        for i in range(len(self.formula.variable)):
            var=self.formula.variable[i]
            ii=self.var.index(var)
            if ii==0:
                self.txt_b_m.insert(0,df['b']['mean']) 
                if self.iProb:
                    self.txt_b_c.insert(0,df['b']['cov'])
            if ii==1:
                self.txt_d_m.insert(0,df['d']['mean'])
                if self.iProb:
                    self.txt_d_c.insert(0,df['d']['cov'])
            if ii==2:
                self.txt_m_m.insert(0,df['m']['mean'])
                if self.iProb:
                    self.txt_m_c.insert(0,df['m']['cov'])
            if ii==3:
                self.txt_Su_m.insert(0,df['Su']['mean'])
                if self.iProb:
                    self.txt_Su_c.insert(0,df['Su']['cov'])
            if ii==4:
                self.txt_Sy_m.insert(0,df['Sy']['mean'])
                if self.iProb:
                    self.txt_Sy_c.insert(0,df['Sy']['cov'])
            if ii==5:
                self.txt_Lsh_m.insert(0,df['Lsh']['mean'])
                if self.iProb:
                    self.txt_Lsh_c.insert(0,df['Lsh']['cov'])
            if ii==6:
                self.txt_th_m.insert(0,df['th']['mean'])
                if self.iProb:
                    self.txt_th_c.insert(0,df['th']['cov'])
            if ii==7:
                self.txt_Limp_m.insert(0,df['Limp']['mean'])
                if self.iProb:
                    self.txt_Limp_c.insert(0,df['Limp']['cov'])
            if ii==8:
                if self.iProb:
                    self.txt_v_m.insert(0,df['v']['mean'])
                    if self.iProb:
                        self.txt_v_c.insert(0,df['v']['cov'])
                else:
                    self.txt_v_m.insert(0,'0')
        for i in range(len(self.formula.const)):
            const=self.formula.const[i]
            ii=self.var.index(const)
            if ii==0:
                self.txt_b_m.insert(0,df['b']['mean'])
            if ii==1:
                self.txt_d_m.insert(0,df['d']['mean'])
            if ii==2:
                self.txt_m_m.insert(0,df['m']['mean'])
            if ii==3:
                self.txt_Su_m.insert(0,df['Su']['mean'])
            if ii==4:
                self.txt_Sy_m.insert(0,df['Sy']['mean'])
            if ii==5:
                self.txt_Lsh_m.insert(0,df['Lsh']['mean'])
            if ii==6:
                self.txt_th_m.insert(0,df['th']['mean'])
            if ii==7:
                self.txt_Limp_m.insert(0,df['Limp']['mean'])
            if ii==8:
                self.txt_v_m.insert(0,df['v']['mean'])
            if ii==9:
                self.txt_ro_imp.insert(0,df['ro_imp']['mean'])
            #if ii==10:
                #self.txt_a10.insert(0,df['a10']['mean'])
            if ii==10:
                self.combo_shape.set(df['shape'])
            if ii==11:
                self.combo_frag.set(df['fragment'])
            if ii==12:
                self.combo_mat.set(df['Material'])

    def dfAdd(self,df):
        #['b', 'd', 'm', 'Su', 'Sy','Lsh', 'th', 'Limp', 'v']
        for i in range(len(self.formula.variable)):
            var=self.formula.variable[i]
            ii=self.var.index(var)
            if ii==0: 
                df['b']={}
                df['b']['mean']=float(self.txt_b_m.get())
                if self.iProb:
                    df['b']['cov']=float(self.txt_b_c.get())
                    df['b']['dist']='normal'  #正規分布と仮定
            if ii==1:
                df['d']={} 
                df['d']['mean']=float(self.txt_d_m.get())
                if self.iProb:
                    df['d']['cov']=float(self.txt_d_c.get())
                    df['d']['dist']='normal'  #正規分布と仮定
            if ii==2:
                df['m']={} 
                df['m']['mean']=float(self.txt_m_m.get())
                if self.iProb:
                    df['m']['cov']=float(self.txt_m_c.get())
                    df['m']['dist']='normal'  #正規分布と仮定
            if ii==3:
                df['Su']={} 
                df['Su']['mean']=float(self.txt_Su_m.get())
                if self.iProb:
                    df['Su']['cov']=float(self.txt_Su_c.get())
                    df['Su']['dist']='normal'  #正規分布と仮定
            if ii==4:
                df['Sy']={} 
                df['Sy']['mean']=float(self.txt_Sy_m.get())
                if self.iProb:
                    df['Sy']['cov']=float(self.txt_Sy_c.get())
                    df['Sy']['dist']='normal'  #正規分布と仮定
            if ii==5:
                df['Lsh']={} 
                df['Lsh']['mean']=float(self.txt_Lsh_m.get())
                if self.iProb:
                    df['Lsh']['cov']=float(self.txt_Lsh_c.get())
                    df['Lsh']['dist']='normal'  #正規分布と仮定
            if ii==6:
                df['th']={} 
                df['th']['mean']=float(self.txt_th_m.get())
                if self.iProb:
                    df['th']['cov']=float(self.txt_th_c.get())
                    df['th']['dist']='normal'  #正規分布と仮定
            if ii==7:
                df['Limp']={} 
                df['Limp']['mean']=float(self.txt_Limp_m.get())
                if self.iProb:
                    df['Limp']['cov']=float(self.txt_Limp_c.get())
                    df['Limp']['dist']='normal'  #正規分布と仮定
            if ii==8:
                df['v']={} 
                df['v']['mean']=float(self.txt_v_m.get())
                if self.iProb:
                    df['v']['cov']=float(self.txt_v_c.get())
                    df['v']['dist']='normal'  #正規分布と仮定

        return df
    def dfAddC(self,df):
        #['b', 'd', 'm', 'Su', 'Sy','Lsh', 'th', 'Limp', 'v']
        for i in range(len(self.formula.const)):
            const=self.formula.const[i]
            ii=self.var.index(const)
            if ii==0:
                df['b']={} 
                df['b']['mean']=float(self.txt_b_m.get())
            if ii==1:
                df['d']={} 
                df['d']['mean']=float(self.txt_d_m.get())
            if ii==2:
                df['m']={} 
                df['m']['mean']=float(self.txt_m_m.get())
            if ii==3:
                df['Su']={} 
                df['Su']['mean']=float(self.txt_Su_m.get())
            if ii==4:
                df['Sy']={} 
                df['Sy']['mean']=float(self.txt_Sy_m.get())
            if ii==5:
                df['Lsh']={} 
                df['Lsh']['mean']=float(self.txt_Lsh_m.get())
            if ii==6:
                df['th']={} 
                df['th']['mean']=float(self.txt_th_m.get())
            if ii==7:
                df['Limp']={} 
                df['Limp']['mean']=float(self.txt_Limp_m.get())
            if ii==8:
                df['v']={} 
                df['v']['mean']=float(self.txt_v_m.get())
            if ii==9:
                df['ro_imp']={} 
                df['ro_imp']['mean']=float(self.txt_ro_imp.get())
            #if ii==10:
                #df['a10']={} 
                #df['a10']['mean']=float(self.txt_a10.get())
            if ii==10:
                df['shape']={}
                df['shape']=self.combo_shape.get()
            if ii==11:
                df['fragment']={}
                df['fragment']=self.combo_frag.get()
            if ii==12:
                df['Material']={}
                df['Material']=self.combo_mat.get()
        return df
    def makeNormal(self):
        self.Clear()
        var=self.formula.variable
        for i in range(len(var)):
            ii=self.var.index(var[i])
            if ii==0:
                self.txt_b_m['state']='normal'
                if self.iProb==1: self.txt_b_c['state']='normal'
            if ii==1:
                self.txt_d_m['state']='normal' 
                if self.iProb==1: self.txt_d_c['state']='normal'
            if ii==2:
                self.txt_m_m['state']='normal' 
                if self.iProb==1: self.txt_m_c['state']='normal'
            if ii==3:
                self.txt_Su_m['state']='normal' 
                if self.iProb==1: self.txt_Su_c['state']='normal'
            if ii==4:
                self.txt_Sy_m['state']='normal' 
                if self.iProb==1: self.txt_Sy_c['state']='normal'
            if ii==5:
                self.txt_Lsh_m['state']='normal' 
                if self.iProb==1: self.txt_Lsh_c['state']='normal'
            if ii==6:
                self.txt_th_m['state']='normal' 
                if self.iProb==1: self.txt_th_c['state']='normal'
            if ii==7:
                self.txt_Limp_m['state']='normal' 
                if self.iProb==1: self.txt_Limp_c['state']='normal'
            if ii==8:
                if self.iProb:
                    self.txt_v_m['state']='normal'                
                    if self.iProb==1: 
                        self.txt_v_c['state']='normal'                    
            if ii==9:
                self.txt_ro_imp['state']='normal'
            #if ii==10:
                #self.txt_a10['state']='normal'
            if ii==10:
                self.combo_shape['state']='normal'
            if ii==11:
                self.combo_frag['state']='normal'
            if ii==12:
                self.combo_mat['state']='normal'
                
            

    def makeConst(self):
        const=self.formula.const #適用範囲チェック用の文字列
        self.lbl_b['text']=''
        self.lbl_d['text']=''
        self.lbl_m['text']=''
        self.lbl_Su['text']=''
        self.lbl_Sy['text']=''
        self.lbl_Lsh['text']=''
        self.lbl_th['text']=''
        self.lbl_Limp['text']=''
        self.lbl_v['text']=''
        self.lbl_ro_imp['text']=''
        #self.lbl_a10['text']=''
        self.lbl_shape['text']=''
        self.lbl_frag['text']=''
        self.lbl_Material['text']=''
        for i in range(len(const)):
            self.setAst(self.var.index(const[i]))       
    def setAst(self,i):
        if i==0: self.lbl_b['text']='*'; self.txt_b_m['state']='normal'
        if i==1: self.lbl_d['text']='*'; self.txt_d_m['state']='normal'
        if i==2: self.lbl_m['text']='*'; self.txt_m_m['state']='normal'
        if i==3: self.lbl_Su['text']='*'; self.txt_Su_m['state']='normal'
        if i==4: self.lbl_Sy['text']='*'; self.txt_Sy_m['state']='normal'
        if i==5: self.lbl_Lsh['text']='*'; self.txt_Lsh_m['state']='normal'
        if i==6: self.lbl_th['text']='*'; self.txt_th_m['state']='normal'
        if i==7: self.lbl_Limp['text']='*'; self.txt_Limp_m['state']='normal'
        if i==8: self.lbl_v['text']='*'; self.txt_v_m['state']='normal'
        if i==9: self.lbl_ro_imp['text']='*' ; self.txt_ro_imp['state']='normal'
        #if i==10: self.lbl_a10['text']='*' ; self.txt_a10['state']='normal' 
        if i==10: self.lbl_shape['text']='*' ; self.combo_shape.config(state='normal');self.combo_shape['values']=['flat','hemisperical']
        if i==11: self.lbl_frag['text']='*' ; self.combo_frag.config(state='normal');self.combo_frag['values']=['Standard','Alternative']
        if i==12: self.lbl_Material['text']='*' ; self.combo_mat.config(state='normal'); self.combo_mat['values']=self.formula.MatList()      
    def makeDisable(self):
        #['b', 'd', 'm', 'Su', 'Sy','Lsh', 'th', 'Limp', 'v']
        self.txt_b_m['state']='disable'
        self.txt_b_c['state']='disable'
        self.txt_d_m['state']='disable'
        self.txt_d_c['state']='disable'
        self.txt_m_m['state']='disable'
        self.txt_m_c['state']='disable'
        self.txt_Su_m['state']='disable'
        self.txt_Su_c['state']='disable'
        self.txt_Sy_m['state']='disable'
        self.txt_Sy_c['state']='disable'
        self.txt_Lsh_m['state']='disable'
        self.txt_Lsh_c['state']='disable'
        self.txt_th_m['state']='disable'
        self.txt_th_c['state']='disable'
        self.txt_Limp_m['state']='disable'
        self.txt_Limp_c['state']='disable'
        self.txt_v_m.insert(0,'0')
        self.txt_v_m['state']='disable'
        self.txt_v_c['state']='disable'
        self.txt_ro_imp.insert(0,'10000')
        #self.txt_a10.insert(0,'1750(aluminum),4000(RHA)')
        self.txt_ro_imp['state']='disable'
        #self.txt_a10['state']='disable'
        #self.txt_shape.insert(0,'flat or hemispherical')
        #self.txt_shape['state']='disable'
        self.combo_shape.config(state='disabled')
        #self.txt_frag.insert(0,'Standard or Alternative')
        #self.txt_frag['state']='disable'
        self.combo_frag.config(state='disabled')
        self.combo_mat.config(state='disabled')
                          
    def change_selected_item(self):
        self.Clear()
        iFormula=self.var_item.get()
        self.iProb=self.var_prob.get()
        if iFormula==0:
            self.formula=pen.AlyLi_M()
        if iFormula==1:
            self.formula=pen.BRL_M()
        if iFormula==2:
            self.formula=pen.DeMarre_M()
        if iFormula==3:
            self.formula=pen.Jowett_M()
        if iFormula==4:
            self.formula=pen.Lambert_M()
        if iFormula==5:
            self.formula=pen.Neilson_M()
        if iFormula==6:
            self.formula=pen.Ohte_M()
        if iFormula==7:
            self.formula=pen.SRI_M()
        if iFormula==8:
            self.formula=pen.SwRI_M()
        if iFormula==9:
            self.formula=pen.THOR_M()
        if iFormula==10:
            self.formula=pen.WenJones_M()
        self.makeDisable()
        self.makeNormal()
        self.makeConst()
        self.txt_Text.delete("1.0","end")
        self.txt_Text.insert(1.0,self.formula.__doc__)
        return
class InterPenet():
    def __init__(self):
        root = tk.Tk()
        root.geometry('1000x700')
        app = Application(master=root)
        app.mainloop()        

