# -*- coding: utf-8 -*-
"""
Created on Sun May 21 21:21:30 2017

@author: CwPc
"""
from scipy import odr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import math
import pickle

from matplotlib import patches as mpatches
from matplotlib import lines as lines





#%%
import pathlength as pl

plotting = False

#%%

class c():
    def __init__(self,name='NoName',**kwargs):#l=None,D=None,fitfunc=None,guess=None,xlabel='x',sigma=None
#        [self.A,self.C],fit_cov = curve_fit(fitfunc,l,D1,guess)
#        print(self.A)
#        print(self.C)
        self.name = name
        self.v = kwargs

        if 'l' in self.v:
            linear = odr.Model(self.v['fitfunc'])
            mydata = odr.RealData(self.v['l'], self.v['D'],sy=self.v['sigma'])
            myodr = odr.ODR(mydata, linear, beta0 = self.v['guess'],maxit=5000)
    #        myodr.set_job(fit_type=2)
            myoutput = myodr.run()
            self.v['A']   = myoutput.beta
            self.A = myoutput.beta
            self.v['fit_cov'] = myoutput.cov_beta
            self.v['fit_sd'] = myoutput.beta
            myoutput.pprint()
            
            print(f'name:{self.name}')
            print(f'A:{self.v["A"]}')
#            self.save_fit()
#        else:
#            self.load_fit()
                

    def show_fit(self):
#        labels = []
#        kvp = [80,100,120,140]
#        for k in kvp:
#            for c in lengths:
#                labels.append(f'{str(c)} mm/{str(k)} kVp')

        lin=np.linspace(0,120,120)
        for i in np.arange(self.D.shape[0]):
            curve, = plt.plot(lin,self.v['fitfunc'](self.v['A'][[0,i+1]],lin),label = labels[i],zorder=1)
            col = curve.get_color()
            plt.errorbar(self.v['l'][i,:],self.v['D'][i,:],color = col, marker ='o',linestyle = 'none',markersize = 5,capsize=2,zorder=2)
        maxval = math.ceil(np.max(self.v['D'])*1.1*10)/10
        plt.axis((0,140,0,maxval))
        plt.ylabel(r'Dose/CTDI$_{vol}$')
        plt.xlabel(self.v['xlabel'])
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

        plt.savefig('out/'+self.name+'.eps',format='eps',dpi=600,bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.show()

        #plt.show()
    def show_discrepancy(self):
        res = self.v['fitfunc'](self.v['A'],self.v['l'])
        res = res/self.v['D']
        self.v['res']=res
#        print(res.std())
#        plt.ylabel(r'$k_{length}$')
        plt.boxplot(res, 0, 'rs',0)
        plt.xlabel('Fit error')
        plt.savefig('out/'+self.name+'_error.eps',format='eps',dpi=600)
        plt.show()

        
    def save_fit(self):
        pass
#        pickle.dump(self.v,open('fit/'+self.name,'wb'))
    
#    def load_fit(self):
#        self.v = pickle.load(open('fit/'+self.name,'rb'))
#        self.A = self.v['A']


class device_settings():
    def __init__(self,dd_df='default',bt_df='default'):
#        This list should come by grace of the input into the app, or manually here
#        bt_list = [(filter_type,kv,fn)]
#        dd_list = [(kv,fn)]
#           
        
        if dd_df == 'default':
            self.dd_list = [[80,'dat/dd/100.xlsx'],
                            [100,'dat/dd/100.xlsx'],
                            [120,'dat/dd/120.xlsx'],
                            [140,'dat/dd/140.xlsx']]
        else:
            for row in dd_df.iterrows():
                print(row)
                
            
            
            
            
            self.dd_list = dd_df
        if bt_df == 'default':
            defaulty = []
            kvs = [80,100,120,140]
            types = ['body','head']
            bt_dir = 'dat/bowtie/'
            for kv in kvs:
                for filt in types:
                    defaulty.append([filt,kv,bt_dir+filt+'.csv'])
            
            self.bt_list = defaulty
        else:
            self.bt_list = bt_df









    




        

#%%
#Fit functions

def fklength(A,X):
    l = X[1]
    kvp = X[2]
    ma=A[1]
    mb = A[2]
#    mc=A[5]
    a=(ma+mb*(kvp-80)/60)
#    t=-mc*(kvp-80)/60
#    return np.log(l*(ma+mb*(kvp-80)/60)+1)
#    return 1-np.exp(-l**.5*a)
#    return 2 - np.exp(-l*(ma+ma*t)) - np.exp(-l*(mb+mb*t))
    return 1-np.exp(-l**.5*a)

#def fkasym(A,X):
#    #X:l,kvp,asym
#    asa,asb=A[0],A[1]
#    l,kvp,asym = X[0],X[1],X[2]
#    return 1/np.cosh((asa+asb*kvp)*l*asym)
    
def fkcoll(A,X):
    #X:[lcoll,l,kvp]
    Aa,Ab=A[3],A[4]
    lcoll,l,kvp=X[0],X[1],X[2]
    return 1/(1+np.log(l/lcoll)/(Aa+Ab*(kvp-80)/60))

#Combined fit function
#Three parameter fit, with C as residual

def Afit(A,X):
    #X:[lcoll,l,kvp]
    #A:[C,lcolla,lcollb,la,lb]
    C=A[0]
    kcoll = fkcoll(A,X)
    klength = fklength(A,X)
    return C * kcoll * klength

#Function for running the fit, returns a fitted calibration class object
def fit_data(df,name):
    Amask = (df.LR==16) & (df.angle==0) & (df.name!='ctdi1402') & (df.type=='ax')
    dfA = df[Amask]
    Aguess = [ 0.858595 ,   0.02192678, 0.00661411,  6.28199231 , 2.68660476]
    cA = c(name=name,l=(dfA.lcoll,dfA.lscan,dfA.kvp),D=dfA.Dnorm,sigma=dfA.Dnormerr, fitfunc = Afit,guess = Aguess, xlabel = 'Scan length (mm)')
    return cA











#%%
#Create new calibration based on measured data
def make_fit(skin_data_fn,name,dd_list='default',bt_list='default'):
    #Load the data set
    df = pd.read_excel(skin_data_fn)

    #create Dnorm column, and corresponding error column
    df['Dnorm'] = df.D/df.ctdi
    df['Dnormerr'] = df.Derr*df.D/df.ctdi

    #Apply the fit function
    cA = fit_data(df,name)
    cDevices = device_settings(dd_list = dd_list,bt_list = bt_list)
    
    
    
    df = apply_fit(df,cA,cDevices)
    
    #If we're calibrating, we can get one final value:
    #Calculate ratio of measured to predicteded
    df['predicted'] = df.D/df.predicted
    pickle.dump([cA,cDevices],open('fit/'+name+'.fit','wb'))
    
    
    plot_fit(df,cA)
    return df,cA,cDevices




def load_fit(fn):
    try:
        
        fit = pickle.load(open(fn,'rb'))
        
    except:
        print('Could not load fit')
    return fit
    
    

#Predict dose based on values in a pandas structure
def apply_fit(df,cA,cDevices):

    #Set the bowtie and dd data in pl
    pl.Devices.set_bt_list(cDevices.bt_list)
    pl.Devices.set_dd_list(cDevices.dd_list)
    
    #Get geometry data from pathlength file:
    df['ksize'] = get_ksize(df)

    df['kangle'] = get_kangle(df)

    
    #Extend length for helical scans to compensate for over-ranging
    m=df.type=='hel'
    t = df.lscan+df.lcoll*df.pitch/2
    t[~m]=df.lscan[~m]
    X=[df.lcoll,t,df.kvp]
    
    #Find length and collimation functions    
    X=[df.lcoll,df.lscan,df.kvp]
    df['klength'] = fklength(cA.A,X)
    df['kcoll'] = fkcoll(cA.A,X)
    
    df['kall'] = cA.A[0]*df.klength*df.kcoll*df.ksize*df.kangle*df.pitch
    
    #Calculate dose prediction
    df['predicted'] = df.ctdi*df.kall

    return df



def get_ksize(idf):
    ksize = idf.apply(lambda x: pl.relative_dose(x['AP'], x['LR'], x['kvp'], x['filter']), axis=1)
    return ksize

def get_kangle(idf):
    kangle = idf.apply(lambda x: pl.relative_dose_to_front(x['AP'], x['LR'], x['kvp'], x['filter'], x['angle']/180*np.pi), axis=1)
    return kangle



def predict(df,fit_fn=None, cA=None,cDevices = None):
    if fit_fn:
        cA,cDevices = load_fit(fit_fn)
    df = apply_fit(df,cA,cDevices)
    return df

def process_spreadsheet(fn,fit_fn = 'fit/ge_optima_660.fit',cA = None, cDevices = None):
    df = pd.read_excel(fn)
    df = predict(df,fit_fn=fit_fn,cA=cA,cDevices=cDevices)
    df.to_excel(fn)


#Only do the rest of this stuff if plotting is enabled
#%%
#plotting=True

def plot_fit(df,cA):
    
    
    #Plot the ratio of predicted to measured values
    masky = (df.angle==0)
    print(df.predicted[masky].std()/df.predicted[masky].mean())
    df = df[masky].sort_values(['LR','kvp','lscan','lcoll']).reset_index()
    sep_gaps = df.LR.value_counts().sort_index().cumsum().reset_index()
    
    seps = [df.index.min()-2]
    for s in sep_gaps.LR[:-1]:
        seps.append(s-0.5)
    
    seps.append(sep_gaps.LR.iloc[-1]+2)
    seps = np.array(seps)
    print(seps)
#    = np.array([-5,38.5,46.5,49.5,55])
    #color to kvp
    #scan length to shape
    kvps = (80,100,120,140)
    colors = ('C0','C1','C2','C3')
    #kvpcol = (1,2,3,4)
    ls = ((20,'x'),(40,'o'),(120,'+'))
    #lshape = ('x','o','+')
    
    fig,ax = plt.subplots()
    
    ax.axhline(y=1, xmin=0, xmax=1, linestyle = '--',color = 'k',zorder=6)
    
    
    xticks = (seps[:-1]+seps[1:])/2
    
    rs = ['16 (Head)','20','32 (CTDI body)','40']
    for i,x in enumerate(seps):
        if i%2==1:
            ax.axvspan(seps[i-1],seps[i], alpha=0.15, color='grey')
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(rs)
    
    ax.set_xlabel('Phantom diameter/type (cm)')
    
    

    for i,kvp in enumerate(kvps):
        for l,shape in ls:
            m = (df.kvp==kvp)&(df.lscan==l) & (df.angle==0)& (df.type=='ax')
            ax.errorbar(df.index[m],df.D[m]/df.predicted[m],df.Derr[m],marker=shape,linestyle='',color=colors[i],zorder =1)
    

    
    testh = [mpatches.Patch(color = colors[i],label = str(kvps[i])+' kVp') for i ,k in enumerate(kvps)]
    testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm scan length') for i ,l in enumerate(ls)]
    
    
    ax.set_xbound(lower=df.index.min()-2,upper = df.index.max()+2)
    ax.set_ybound(lower = .7,upper = 1.3)
    ax.legend(handles = testh+testi,ncol=2,fancybox=True)
    #ax.set_ylim([.65,1.35])
    ax.set_ylabel('Measured/predicted surface dose')
    
    fig.tight_layout()
    fig.savefig('out/'+'datavsmodel.png',format='png',dpi=300)
    fig.savefig('out/'+'datavsmodel.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    
def plot_pre_fit(df,cA):
    
    
    #Plot the ratio of predicted to measured values
    masky = (df.angle==0)
    print(df.predicted[masky].std()/df.predicted[masky].mean())
    df = df[masky].sort_values(['LR','kvp','lscan','lcoll']).reset_index()
    sep_gaps = df.LR.value_counts().sort_index().cumsum().reset_index()
    
    seps = [df.index.min()-2]
    for s in sep_gaps.LR[:-1]:
        seps.append(s-0.5)
    
    seps.append(sep_gaps.LR.iloc[-1]+2)
    seps = np.array(seps)
    print(seps)
#    = np.array([-5,38.5,46.5,49.5,55])
    #color to kvp
    #scan length to shape
    kvps = (80,100,120,140)
    colors = ('C0','C1','C2','C3')
    #kvpcol = (1,2,3,4)
    ls = ((20,'x'),(40,'o'),(120,'+'))
    #lshape = ('x','o','+')
    
    fig,ax = plt.subplots()
    
#    ax.axhline(y=1, xmin=0, xmax=1, linestyle = '--',color = 'k',zorder=6)
    
    
    xticks = (seps[:-1]+seps[1:])/2
    
    rs = ['16 (Head)','20','32 (CTDI body)','40']
    for i,x in enumerate(seps):
        if i%2==1:
            ax.axvspan(seps[i-1],seps[i], alpha=0.15, color='grey')
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(rs)
    
    ax.set_xlabel('Phantom diameter/type (cm)')
    
    

    for i,kvp in enumerate(kvps):
        for l,shape in ls:
            m = (df.kvp==kvp)&(df.lscan==l) & (df.angle==0)& (df.type=='ax')
            ax.errorbar(df.index[m],df.D[m]/df.ctdi[m],df.Derr[m],marker=shape,linestyle='',color=colors[i],zorder =1)
    

    
    testh = [mpatches.Patch(color = colors[i],label = str(kvps[i])+' kVp') for i ,k in enumerate(kvps)]
    testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm scan length') for i ,l in enumerate(ls)]
    
    
    ax.set_xbound(lower=df.index.min()-2,upper = df.index.max()+2)
    ax.set_ybound(lower = .7)
    ax.legend(handles = testh+testi,ncol=2,fancybox=True)
    #ax.set_ylim([.65,1.35])
    ax.set_ylabel('Measured dose/CTDIvol')
    
    fig.tight_layout()
#    fig.savefig('out/'+'datavsmodelpre.png',format='png',dpi=300)
    fig.savefig('out/'+'datavsmodelpre.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    
    
    
    #print((df.predicted[Amask]).std()/(df.predicted[Amask]).mean())
    
#a = make_fit('dat/ge_optima_660_raw.xlsx','ge_optima_660')
#    return df
    
#plot_fit(df,cA)

    #%%

def plotting(df,cA):
    df['Dnorm'] = df.D/df.ctdi
    #Plot the raw data 
    
    fig,ax = plt.subplots()
    
    kvps = (80,100,120,140)
    
    colls = (20,40,120)
    colors = ('C0','C1','C2','C3')
    #kvpcol = (1,2,3,4)
    ls = ((10,'x'),(20,'o'),(40,'+'))
    
    
    
    
    seps = np.array([-5,7.5,15.5,23.5,36])
    xticks = (seps[:-1]+seps[1:])/2
    
    rs = ['80 kVp','100 kVp','120 kVp','140 kVp']
    for i,x in enumerate(seps):
        if i%2==1:
            ax.axvspan(seps[i-1],seps[i], alpha=0.15, color='grey')
    
    
    ax.set_xlim((-2,36))
    ax.set_ylim((.9,2))
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(rs)
    ax.set_xlabel('Beam energy')
    
    
    
    
    
    masky = (df.LR==16) & (df.angle==0)& (df.type=='ax')
    testdf = df[masky].sort_values(by=['kvp','lscan','lcoll'])
    
    for i,coll in enumerate(colls):
        for l,shape in ls:
            m = (testdf.lscan==coll)&(testdf.lcoll==l)
            ax.errorbar(testdf.reset_index().index[m],testdf.Dnorm[m],testdf.Dnorm[m]*testdf.Derr[m],marker=shape,linestyle='',color=colors[i],zorder =1)
    
    #        plt.plot()
    
    testh = [mpatches.Patch(color = colors[i],label = str(colls[i])+' mm scan length') for i ,k in enumerate(colls)]
    testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm coll.') for i ,l in enumerate(ls)]
    
    ax.legend(handles = testh+testi,ncol=2,fancybox=True)
    #ax.set_ylim([.65,1.3])
    ax.set_ylabel(r'Surface dose/CTDI$_{vol}$')
    
    fig.tight_layout()
    fig.savefig('out/'+'ctdidata.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    
    

    
    #Plot klength, alonside a few data points
    masky = (df.angle==0)&(df.LR==16) & (df.name!='ctdi1402')& (df.type=='ax')#&(df.l==df.lcoll)
    
    x = np.linspace(0,120,120)
    
    fig,ax = plt.subplots()
    
    kvs = [80,100,120,140]
    colors = [0,1,2,3]
    for i, kv in enumerate(kvs):
        y = fklength((cA.A),(x,x,kv))
        ax.plot(x,y,label = str(kv)+' kVp',color='C'+str(i))
        m=masky&(df.kvp==kv)
        ax.errorbar(df.lscan[m],df.Dnorm[m]/df.kcoll[m]/cA.A[0],(df.Derr*df.Dnorm/cA.A[0])[m],fmt='o',color='C'+str(i),alpha=.6,label='_nolegend_')
        
    ax.legend()
    ax.set_ylabel(r'Scan length correction factor $k_{length}$')
    ax.set_xlabel('Scan length (mm)')
    fig.tight_layout()
    fig.savefig('out/'+'klength.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    

    
    
    
    
    

    #Plot kcoll, it doesn't really make sense to add data points to this plot
    
    masky = (df.angle==0)&(df.lscan==df.lcoll)&(df.LR==16)
    
    x = np.linspace(1,16,50)
    
    fig,ax = plt.subplots()
    
    kvs = [80,100,120,140]
    for i, kv in enumerate(kvs):
        y = fkcoll((cA.A),(1,x,kv))
        ax.plot(x,y,label = str(kv)+' kVp',color='C'+str(i))
        m=masky&(df.kvp==kv)
    #    ax.errorbar(df.l[m],df.Dnorm[m]/cA.A[0],(df.Derr*df.Dnorm/cA.A[0])[m],fmt='o',color='C'+str(i),alpha=.6)
    
    #ax.set_ylim((0,1))
    ax.legend()
    ax.set_ylabel(r'Collimation correction factor $k_{coll}$')
    ax.set_xlabel('Number of rotations (n)')
    fig.tight_layout()
    fig.savefig('out/'+'kcoll.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    
    

    #Plot helical mode points
    
    #plt.plot(df.predicted[m]/cA.A[0],'o')
    
    #Calculated over-rang for helical scans
    m=df.type=='hel'
    
#    X=[df.lcoll,t,df.kvp]
#    df['klength'] = fklength(cA.A,X)
    
    
    
    #color to kvp
    #scan length to shape
    kvps = (80,100,120,140)
    es = (0,50,75)
    colors = ('C0','C1','C2')
    #kvpcol = (1,2,3,4)
    ls = ((20,'x'),(40,'o'))
    #lshape = ('x','o','+')
    
    fig,ax = plt.subplots()
    
    ax.axhline(y=1, xmin=0, xmax=1, linestyle = '--',color = 'k',zorder=6)
    #
    seps = np.array([53.5,55.5,57.5,59.5,61.5])
    xticks = (seps[:-1]+seps[1:])/2
    #
    rs = ['80 kVp','100 kVp','120 kVp','140 kVp']
    for i,x in enumerate(seps):
        if i%2==1:
            ax.axvspan(seps[i-1],seps[i], alpha=0.1, color='grey')
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(rs)
    
    ax.set_xlabel('Peak beam energy (kV)')
    
    
    ax.set_xlim((53.5,61.5))
    
    for i,e in enumerate(es):
        t = df.lscan+df.lcoll*df.pitch*e/100
        X = [df.lcoll,t,df.kvp]
        klength = fklength(cA.A,X)
        predicted=df.Dnorm/klength/df.kcoll/df.ksize/df.pitch
        
        for l,shape in ls:
            m = (df.lcoll==l) & (df.angle==0)& (df.type=='hel')
            ax.errorbar(df.index[m],predicted[m]/cA.A[0],df.Derr[m]*0,marker=shape,linestyle='',color='C'+str(i),zorder =1)
    
    #        plt.plot()
    
    testh = [mpatches.Patch(color = colors[i],label = 'e = '+str(es[i])+'%') for i ,k in enumerate(colls)]
    testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm collimation') for i ,l in enumerate(ls)]
    
    
    
    ax.legend(handles = testh+testi,ncol=2,fancybox=True)
    ax.set_ylim([.65,1.35])
    ax.set_ylabel('Measured/predicted surface dose')
    
    fig.tight_layout()
    fig.savefig('out/'+'hel.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    
    
    
    
    
    #%%
    #Check against the delas heras data set
    
    delas = pd.read_excel('delas.xlsx')
    delas['ksize'] = delas['LR']/delas['LR']
    #delas['ksize'] = delas.apply(lambda x: pl.relative_dose(x['LR'], x['AP'], x['kvp'], x['filt']), axis=1)
    delas['predicted'] = Afit(cA.A,[delas.lcoll,delas.lscan,delas.kvp])
    delas['ratio'] = delas.D/delas.ctdi/delas.predicted
    delas.to_excel('delas.xlsx')
    plt.plot()
    
    #%%
    







#%%
default_fn = 'fit/ge_optima_660.fit'
cA,cDevices = load_fit(default_fn)
df = pd.read_excel('dat/ge_optima_660_raw.xlsx')
process_spreadsheet('predict.xlsx',cA=cA,cDevices=cDevices)

#df = predicted(df,cA=cA,cDevices=cDevices)
#
#plot_fit(df,cA)
#plotting(df,cA)








