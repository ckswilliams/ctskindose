# -*- coding: utf-8 -*-
"""
Created on Sun May 21 21:21:30 2017

@author: CwPc
"""
from scipy import odr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import pickle
import itertools

from matplotlib import patches as mpatches
from matplotlib import lines as lines

import pathlength as pl

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
            myodr.set_job(fit_type=2)
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
                

    # def show_fit(self):
    #     labels = []
    #     kvp = [80,100,120,140]
    #     for k in kvp:
    #         for c in self.v['l']:
    #             labels.append(f'{str(c)} mm/{str(k)} kVp')
    #
    #     lin=np.linspace(0,120,120)
    #     for i in np.arange(self.D.shape[0]):
    #         curve, = plt.plot(lin,self.v['fitfunc'](self.v['A'][[0,i+1]],lin),label = labels[i],zorder=1)
    #         col = curve.get_color()
    #         plt.errorbar(self.v['l'][i,:],self.v['D'][i,:],color = col, marker ='o',linestyle = 'none',markersize = 5,capsize=2,zorder=2)
    #     maxval = math.ceil(np.max(self.v['D'])*1.1*10)/10
    #     plt.axis((0,140,0,maxval))
    #     plt.ylabel(r'Dose/CTDI$_{vol}$')
    #     plt.xlabel(self.v['xlabel'])
    #     lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    #
    #     plt.savefig('out/'+self.name+'.eps',format='eps',dpi=600,bbox_extra_artists=(lgd,), bbox_inches='tight')
    #     plt.show()

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


class DeviceSettings():
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



from scipy.optimize import curve_fit
class CF():
    def __init__(self, func, D, X, guess, sigma=None):
        self.A, self.v = curve_fit(func, X, D, guess, sigma, absolute_sigma=True)
        

#%% Original fit functions
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
    return 1-np.exp(-l*a)

#def fkasym(A,X):
#    #X:l,kvp,asym
#    asa,asb=A[0],A[1]
#    l,kvp,asym = X[0],X[1],X[2]
#    return 1/np.cosh((asa+asb*kvp)*l*asym)
    
def fkcoll(A,X):
    #X:[lcoll,l,kvp]
    Aa,Ab=A[3],A[4]
    lcoll,lscan,kvp=X[0],X[1],X[2]
    return 1/(1+np.log(lscan/lcoll)/(Aa+Ab*(kvp-80)/60))




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
    Amask = (df.LR==16) & (df.angle==0) & (str(df.name) !='ctdi1402') & (df.scan_type=='ax')
    dfA = df[Amask]
    Aguess = [ 2.48866957,  0.13324515, -0.01845219,  7.47206539,  0.83062203]
    cA = c(name=name,l=(dfA.lcoll,dfA.lscan,dfA.kvp),D=dfA.Dnorm,sigma=dfA.Dnormerr, fitfunc = Afit,guess = Aguess, xlabel = 'Scan length (mm)')
    return cA


#%% Testing fit functions


def fklength(A,X):
    '''
    kVp-dependent fit function for total scan length
    X:[l0,l,kvp] (beam width, scan length, kVp)
    A:[c, c_penumbra, c_length, c_kvp]
    '''
    l = X[1]
    c_penumbra = A[1]
    c_length = A[2]
    #return 1 + c_length*l
    return 1-0.5*np.exp(-c_length*(l))


#def fkasym(A,X):
#    #X:l,kvp,asym
#    asa,asb=A[0],A[1]
#    l,kvp,asym = X[0],X[1],X[2]
#    return 1/np.cosh((asa+asb*kvp)*l*asym)
    
def fkcoll(A,X):
    '''
    fit function for collimation
    X:[l0,l,kvp] (beam width, scan length, kVp)
    A:[c, c_penumbra, c_length, c_kvp]
    '''
    c_penumbra = A[1]
    lcoll=X[0]
    #return 1/(1+np.log(lscan/lcoll)/(Aa+Ab*(kvp-80)/60))
    return lcoll / (lcoll + c_penumbra)


def fkkvp(A,X):
    '''
    fit function for kvp
    X:[l0,l,kvp] (beam width, scan length, kVp)
    A:[c, c_penumbra, c_length, c_kvp]
    '''
    kvp = X[2]
    c_kvp = A[3]
    return 1-c_kvp*(kvp-70)/140

def Afit(A,*X):
    #X:[lcoll,l,kvp]
    #A:[C,lcolla,lcollb,la,lb]
    C=A[0]
    kcoll = fkcoll(A,X)
    klength = fklength(A,X)
    kkvp = fkkvp(A,X)
    return C * kcoll * klength * kkvp

def fk_curvefit(X,*A):
    #X:[lcoll,l,kvp]
    #A:[C,lcolla,lcollb,la,lb]
    C=A[0]
    kcoll = fkcoll(A,X)
    klength = fklength(A,X)
    kkvp = fkkvp(A,X)
    return C * kcoll * klength * kkvp

#Function for running the fit, returns a fitted calibration class object
def fit_data(df,name):
    Amask = (df.LR==16) & (df.angle==0) & (df.scan_type=='ax')
    dfA = df[Amask]
    #Aguess = [ 2.48866957,  0.13324515, -0.01845219,  7.47206539,  0.83062203]
    Aguess = (1.6,4,.1,.1)# 4 parameter fit instead of 5
    X = (dfA.lcoll,dfA.effective_lscan,dfA.kvp)
    D = dfA.Dnorm
    cA = CF(fk_curvefit, D,X,Aguess, dfA.Dnormerr)
    #cA = c(name=name,l=X,D=dfA.Dnorm,sigma=dfA.Dnormerr, fitfunc = Afit,guess = Aguess, xlabel = 'Scan length (mm)')
    return cA


#%%
#Create new calibration based on measured data
def make_fit(skin_data_fn,name,dd_list='default',bt_list='default'):
    #Load the data set
    
    
    df = pd.read_excel(skin_data_fn)

    #create Dnorm column, and corresponding error column
    df['Dnorm'] = df.D/df.ctdi
    df['Dnormerr'] = df.Derr*df.D/df.ctdi
    df['effective_lscan'] = df.lscan + (df.scan_type=='hel')*df.lcoll*df.pitch/2

    #Apply the fit function
    cA = fit_data(df,name)
    cDevices = DeviceSettings()
    
    
    
    df = apply_fit(df,cA,cDevices)
    
    #If we're calibrating, we can get one final value:
    #Calculate ratio of measured to predicteded

    pickle.dump([cA,cDevices],open('fit/'+name+'.fit','wb'))
    
    
    plot_fit(df)
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
    
    df['Dnorm'] = df.D/df.ctdi
    df['Dnormerr'] = df.Derr*df.D/df.ctdi
    
    #Get geometry data from pathlength file:
    df['ksize'] = get_ksize(df)

    df['kangle'] = get_kangle(df)

    
    #Extend length for helical scans to compensate for over-ranging
    df['effective_lscan'] = df.lscan + (df.scan_type=='hel')*df.lcoll*df.pitch/2

    #Find length and collimation functions    
    X=[df.lcoll,df.effective_lscan,df.kvp]
    df['klength'] = fklength(cA.A,X)
    df['kcoll'] = fkcoll(cA.A,X)
    df['kkvp'] = fkkvp(cA.A,X)
    
    df['kall'] = cA.A[0]*df.klength*df.kcoll*df.ksize*df.kangle*df.pitch*df.kkvp
    
    #Calculate dose prediction
    df['predicted'] = df.ctdi*df.kall
    if 'D' in df.columns:
        df['ratio'] = df.D/df.predicted

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

  

def plot_fit(df):
    
    fig, ax = skin_dose_plot(df, 'ratio')
    
    # Set 
    #ax.set_xbound(lower=df.index.min()-2,upper = df.index.max()+2)
#    ax.set_ybound(lower = .7,upper = 1.3)

    #ax.set_ylim([.65,1.35])
    ax.set_ylabel('Measured/predicted surface dose')
    
    fig.tight_layout()
    #fig.savefig('out/'+'datavsmodel.png',format='png',dpi=300)
    fig.savefig('out/'+'datavsmodel.svg',format='svg',dpi=600)
    plt.show()
    plt.close()
    
    
    
def skin_dose_plot(df, plot_col = 'ratio', error='Dnormerr'):

    #Plot the ratio of predicted to measured values   
    df = df.sort_values(['LR','kvp','lscan','lcoll']).reset_index()
    


    # Get the key values from the dataframe, and assign a color to each kVp
    # Then assign a shape to each scan lenght/beam collimation combo
    kvps = df.kvp.unique()
    kvps.sort()
        
    colors = ('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10')
    shapes = ('*','s','x','o','^','<','>','v','P','1','2','3','4')
    
    # Create iterables with unique scan/beam width with shape
    l_shapes = [(a,b,c) for (a,b),c in zip(df.loc[:,['lcoll','lscan']].drop_duplicates().sort_values(['lcoll','lscan']).values, shapes)]
    
    # Iterable relating kvp to color
    kvp_colors = list(zip(kvps, colors))
    
    # Start the plot
    fig,ax = plt.subplots()
    
    # Add info relating to phantom diameter via the X axis
    ax.axhline(y=1, xmin=0, xmax=1, linestyle = '--',color = 'k',zorder=6)
    

    # Measure where the different phantom size data should is in the index
    sep_gaps = df.LR.value_counts().sort_index().cumsum().reset_index()
    phantom_sizes = {8:'16\n(CTDI head)',16:'32\n(CTDI body)'}
    sep_gaps['gap_name'] = sep_gaps['index'].apply(lambda x: phantom_sizes[x] if x in phantom_sizes else x*2)
    
    seps = [df.index.min()-2]
    for s in sep_gaps.LR[:-1]:
        seps.append(s-0.5)
    
    seps.append(sep_gaps.LR.iloc[-1]+2)
    seps = np.array(seps)
    
    xticks = (seps[:-1]+seps[1:])/2
    
    
    for i,x in enumerate(seps):
        if i%2==1:
            ax.axvspan(seps[i-1],seps[i], alpha=0.15, color='grey')
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(sep_gaps.gap_name)
    
    ax.set_xlabel('Phantom diameter/type (cm)')
    
    for (lcoll, lscan, shape), (kvp, color) in itertools.product(l_shapes, kvp_colors):
        m = (df.kvp==kvp) & (df.lscan==lscan) & (df.lcoll==lcoll)
        if error:
            ax.errorbar(df[m].index, df.loc[m,plot_col], df.loc[m,error],marker=shape, linestyle='', color=color, zorder=1)
        ax.plot(df.index[m], df.loc[m,plot_col], marker=shape,linestyle='',color=color,zorder =1)
    

    # Add color patches and marker meaning to the legend manually
    legend_kvp_colors = [mpatches.Patch(color = color,label = f'{kvp} kVp') for kvp, color in kvp_colors]
    legend_shapes = [
        lines.Line2D([], [],color = 'k',linestyle='', marker=shape,
                     label = f'BW:{lcoll}, SL: {lscan}'
                     ) for lcoll, lscan, shape in l_shapes]
    
    # Set 
    #ax.set_xbound(lower=df.index.min()-2,upper = df.index.max()+2)
#    ax.set_ybound(lower = .7,upper = 1.3)
    ax.legend(handles = legend_kvp_colors+legend_shapes,ncol=2,fancybox=True)
    return fig, ax
    
    
def plot_pre_fit(df):

    fig, ax = skin_dose_plot(df, 'Dnorm')
    #ax.set_ylim([.65,1.35])
    ax.set_ylabel('Measured surface dose/Reported CTDIvol')
    
    fig.tight_layout()
    fig.savefig('out/'+'datavsmodelpre.svg',format='svg',dpi=600)
    #fig.savefig('out/'+'datavsmodel.pdf',format='pdf',dpi=600)
    plt.show()
    plt.close()
    
#plot_fit(df,cA)

    #%% PLOTTING





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
    
    
    
    
    
    masky = (df.LR==16) & (df.angle==0)& (df.scan_type=='ax')
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
    masky = (df.angle==0)&(df.LR==16) & (str(df.name) !='ctdi1402')& (df.scan_type=='ax')#&(df.l==df.lcoll)
    
    x = np.linspace(0,120,120)
    
    fig,ax = plt.subplots()
    
    kvs = [80,100,120,140]
    colors = [0,1,2,3]
    for i, kv in enumerate(kvs):
        y = fklength((cA.A),(x,x,kv))
        ax.plot(x,y,label = str(kv)+' kVp',color='C'+str(i))
        m=masky&(df.kvp==kv)
#        ax.errorbar(df.lscan[m],df.Dnorm[m]/df.kcoll[m]/cA.A[0],(df.Derr*df.Dnorm/cA.A[0])[m],fmt='o',color='C'+str(i),alpha=.6,label='_nolegend_')
        
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
    m=df.scan_type=='hel'
    
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
            m = (df.lcoll==l) & (df.angle==0)& (df.scan_type=='hel')
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
    
    
    
    
    

    #Check against the delas heras data set
    
    delas = pd.read_excel('delas.xlsx')
    delas['ksize'] = delas['LR']/delas['LR']
    #delas['ksize'] = delas.apply(lambda x: pl.relative_dose(x['LR'], x['AP'], x['kvp'], x['filt']), axis=1)
    delas['predicted'] = Afit(cA.A,[delas.lcoll,delas.lscan,delas.kvp])
    delas['ratio'] = delas.D/delas.ctdi/delas.predicted
    delas.to_excel('delas.xlsx')
    plt.plot()
    

#%%
#default_fn = 'fit/ge_optima_660.fit'
#cA,cDevices = load_fit(default_fn)
#df = pd.read_excel('dat/ge_optima_660_raw.xlsx')
#process_spreadsheet('predict.xlsx',cA=cA,cDevices=cDevices)

#df = predicted(df,cA=cA,cDevices=cDevices)
#

#plotting(df,cA)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--fit', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--fit_name', default='fit/combined.fit')
    parser.add_argument('--process', default='dat/combined_data_input.xlsx')
    parser.add_argument('--output_fn', default='out/prediction_output.xlsx')
    args = parser.parse_args()

    if args.fit:
        df, cA, cDev = make_fit('dat/combined_data_input.xlsx','combined')
        cA,cDevices = load_fit('fit/combined.fit')

        df2 = pd.read_excel('dat/combined_data_input.xlsx')
        df2 = apply_fit(df2, cA, cDevices)

        if args.plot:
            plot_fit(df2)
            plot_pre_fit(df2)

    if args.process or not args.fit:
        df = pd.read_excel(args.process)
        cA, cDevices = load_fit(args.fit_name)
        df = apply_fit(df2, cA, cDevices)
        df.to_excel(args.output_fn)






#%%


