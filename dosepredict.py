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

from matplotlib import patches as mpatches
from matplotlib import lines as lines
#%%
import pathlength as pl

#%%

class c():
    def __init__(self,l,D,fitfunc,name='NoName',guess=None,xlabel='x',sigma=None):
#        [self.A,self.C],fit_cov = curve_fit(fitfunc,l,D1,guess)
#        print(self.A)
#        print(self.C)
        self.name = name
        self.l = l
        self.D=D
        self.fitfunc = fitfunc
        self.xlabel = xlabel
        self.sigma = sigma
        linear = odr.Model(self.fitfunc)
        mydata = odr.RealData(l, D,sy=sigma)
        myodr = odr.ODR(mydata, linear, beta0 = guess,maxit=5000)
#        myodr.set_job(fit_type=2)
        myoutput = myodr.run()
        self.A   = myoutput.beta
        self.fit_cov = myoutput.cov_beta
        self.fit_sd = myoutput.beta
        myoutput.pprint()
        
        print(f'name:{self.name}')
        print(f'A:{self.A}')


    def show_fit(self):
#        labels = []
#        kvp = [80,100,120,140]
#        for k in kvp:
#            for c in lengths:
#                labels.append(f'{str(c)} mm/{str(k)} kVp')

        lin=np.linspace(0,120,120)
        for i in np.arange(self.D.shape[0]):
            curve, = plt.plot(lin,self.fitfunc(self.A[[0,i+1]],lin),label = labels[i],zorder=1)
            col = curve.get_color()
            plt.errorbar(self.l[i,:],self.D[i,:],color = col, marker ='o',linestyle = 'none',markersize = 5,capsize=2,zorder=2)
        maxval = math.ceil(np.max(self.D)*1.1*10)/10
        plt.axis((0,140,0,maxval))
        plt.ylabel(r'Dose/CTDI$_{vol}$')
        plt.xlabel(self.xlabel)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

        plt.savefig(self.name+'.eps',format='eps',dpi=600,bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.show()

        #plt.show()
    def show_discrepancy(self):
        res = self.fitfunc(self.A,self.l)
        res = res/self.D
        res=res
        self.res=res
#        print(res.std())
#        plt.ylabel(r'$k_{length}$')
        plt.boxplot(res, 0, 'rs',0)
        plt.xlabel('Fit error')
        plt.savefig(self.name+'_error.eps',format='eps',dpi=600)
        plt.show()

#%%
#Load the data set
df = pd.read_excel('stats.xlsx')
#create Dnorm column, and corresponding error column
df['Dnorm'] = df.D/df.ctdi
df['Dnormerr'] = df.Derr*df.D/df.ctdi


#%%   
#Get geometry data from pathlength file:
df['ksize'] = df.apply(lambda x: pl.relative_dose(x['r'], x['r'], x['kvp'], x['filt']), axis=1)
df['kangle'] = df.apply(lambda x: pl.relative_dose_to_front(x['r'], x['r'], x['kvp'], x['filt'], x['angle']/180*np.pi), axis=1)

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
    #X:[l0,l,kvp]
    Aa,Ab=A[3],A[4]
    l0,l,kvp=X[0],X[1],X[2]
    return 1/(1+np.log(l/l0)/(Aa+Ab*(kvp-80)/60))
#%%
#Fit to the measured data

#We only want to fit to measurements at r=16 and angle =0
Amask = (df.r==16) & (df.angle==0) & (df.name!='ctdi1402') & (df.type=='ax')
dfA = df[Amask]

#Three parameter fit, with C as residual

def Afit(A,X):
    #X:[l0,l,kvp]
    #A:[C,l0a,l0b,la,lb]
    C=A[0]
    kcoll = fkcoll(A,X)
    klength = fklength(A,X)
    return C * kcoll * klength
AA=[.4,.9,-.22,8,-.914]
AA2 = [ 0.858595 ,   0.02192678, 0.00661411,  6.28199231 , 2.68660476]

cA = c((dfA.l0,dfA.l,dfA.kvp),dfA.Dnorm,sigma=dfA.Dnormerr, fitfunc = Afit,guess = AA2, xlabel = 'Scan length (mm)')
#dcA.show_fit()	

X=[df.l0,df.l,df.kvp]
df['klength'] = fklength(cA.A,X)
df['kcoll'] = fkcoll(cA.A,X)


#%%

#Compenaste for helical mode

#Calculated over-rang for helical scans
m=df.type=='hel'
t = df.l+df.l0*df.p/2
t[~m]=df.l[~m]
X=[df.l0,t,df.kvp]
df['klength'] = fklength(cA.A,X)





df['final']=df.Dnorm/df.klength/df.kcoll/df.ksize/df.p

#%%
masky = df.angle==0
print(df.final[masky].std()/df.final[masky].mean())

#color to kvp
#scan length to shape
kvps = (80,100,120,140)
colors = ('C0','C1','C2','C3')
#kvpcol = (1,2,3,4)
ls = ((20,'x'),(40,'o'),(120,'+'))
#lshape = ('x','o','+')

fig,ax = plt.subplots()

ax.axhline(y=1, xmin=0, xmax=1, linestyle = '--',color = 'k',zorder=6)

seps = np.array([-5,38.5,46.5,49.5,55])
xticks = (seps[:-1]+seps[1:])/2

rs = ['32 (CTDI body)','16 (Head)','20','40']
for i,x in enumerate(seps):
    if i%2==1:
        ax.axvspan(seps[i-1],seps[i], alpha=0.15, color='grey')

ax.set_xticks(xticks)
ax.set_xticklabels(rs)

ax.set_xlabel('Phantom diameter/type (cm)')


ax.set_xlim((-2,54))
for i,kvp in enumerate(kvps):
    for l,shape in ls:
        m = (df.kvp==kvp)&(df.l==l) & (df.angle==0)& (df.type=='ax')
        ax.errorbar(df.index[m],df.final[m]/cA.A[0],df.Derr[m],marker=shape,linestyle='',color=colors[i],zorder =1)

#        plt.plot()

testh = [mpatches.Patch(color = colors[i],label = str(kvps[i])+' kVp') for i ,k in enumerate(kvps)]
testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm scan length') for i ,l in enumerate(ls)]

ax.legend(handles = testh+testi,ncol=2,fancybox=True)
#ax.set_ylim([.65,1.35])
ax.set_ylabel('Measured/predicted surface dose')

fig.tight_layout()
fig.savefig('datavsmodel.png',format='png',dpi=300)
fig.savefig('datavsmodel.pdf',format='pdf',dpi=600)
plt.show()
plt.close()
print((df.final[Amask]).std()/(df.final[Amask]).mean())
#%%

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





masky = (df.r==16) & (df.angle==0)& (df.type=='ax')
testdf = df[masky].sort_values(by=['kvp','l','l0'])

for i,coll in enumerate(colls):
    for l,shape in ls:
        m = (testdf.l==coll)&(testdf.l0==l)
        ax.errorbar(testdf.reset_index().index[m],testdf.Dnorm[m],testdf.Dnorm[m]*testdf.Derr[m],marker=shape,linestyle='',color=colors[i],zorder =1)

#        plt.plot()

testh = [mpatches.Patch(color = colors[i],label = str(colls[i])+' mm scan length') for i ,k in enumerate(colls)]
testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm coll.') for i ,l in enumerate(ls)]

ax.legend(handles = testh+testi,ncol=2,fancybox=True)
#ax.set_ylim([.65,1.3])
ax.set_ylabel(r'Surface dose/CTDI$_{vol}$')

fig.tight_layout()
fig.savefig('ctdidata.pdf',format='pdf',dpi=600)
plt.show()
plt.close()
print((df.final[Amask]).std()/(df.final[Amask]).mean())

#%%

masky = (df.angle==0)&(df.r==16) & (df.name!='ctdi1402')& (df.type=='ax')#&(df.l==df.l0)

x = np.linspace(0,120,120)

fig,ax = plt.subplots()

kvs = [80,100,120,140]
colors = [0,1,2,3]
for i, kv in enumerate(kvs):
    y = fklength((cA.A),(x,x,kv))
    ax.plot(x,y,label = str(kv)+' kVp',color='C'+str(i))
    m=masky&(df.kvp==kv)&(df.l==df.l0)
    ax.errorbar(df.l[m],df.Dnorm[m]/df.kcoll[m]/cA.A[0],(df.Derr*df.Dnorm/cA.A[0])[m],fmt='o',color='C'+str(i),alpha=.6,label='_nolegend_')
    
ax.legend()
ax.set_ylabel(r'Scan length correction factor $k_{length}$')
ax.set_xlabel('Scan length (mm)')
fig.savefig('klength.pdf',format='pdf',dpi=600)
plt.show()
plt.close()

#%%





#%%


masky = (df.angle==0)&(df.l==df.l0)&(df.r==16)

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
fig.savefig('kcoll.pdf',format='pdf',dpi=600)
plt.show()
plt.close()


#%%
#Plot helical mode points

#plt.plot(df.final[m]/cA.A[0],'o')

#Calculated over-rang for helical scans
m=df.type=='hel'
X=[df.l0,t,df.kvp]
df['klength'] = fklength(cA.A,X)



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
    t = df.l+df.l0*df.p*e/100
    X = [df.l0,t,df.kvp]
    klength = fklength(cA.A,X)
    final=df.Dnorm/klength/df.kcoll/df.ksize/df.p
    
    for l,shape in ls:
        m = (df.l0==l) & (df.angle==0)& (df.type=='hel')
        ax.errorbar(df.index[m],final[m]/cA.A[0],df.Derr[m]*0,marker=shape,linestyle='',color='C'+str(i),zorder =1)

#        plt.plot()

testh = [mpatches.Patch(color = colors[i],label = 'e = '+str(es[i])+'%') for i ,k in enumerate(colls)]
testi = [lines.Line2D([], [],color = 'k',linestyle='',marker=l[1],label = str(l[0]) + ' mm collimation') for i ,l in enumerate(ls)]



ax.legend(handles = testh+testi,ncol=2,fancybox=True)
ax.set_ylim([.65,1.35])
ax.set_ylabel('Measured/predicted surface dose')

fig.tight_layout()
fig.savefig('hel.pdf',format='pdf',dpi=600)
plt.show()
plt.close()





#%%


delas = pd.read_excel('delas.xlsx')
delas['ksize'] = delas.apply(lambda x: pl.relative_dose(x['r'], x['r'], x['kvp'], x['filt']), axis=1)
delas['final'] = delas.D/delas.ctdi/Afit(cA.A,[delas.l0,delas.l,delas.kvp])
delas['predict'] = Afit(cA.A,[delas.l0,delas.l,delas.kvp])
delas.to_excel('delas.xlsx')

#%%
















#%%
fwfghads
def fit(A,X):
    #X:[l0,l,kvp]
    #A:[C,ma,mb,Aa,Ab]
    C=A[0]
    kcoll = fkcoll(A[1:3],X)
    klength = fklength(A[3:])
    return C * kcoll * klength

    


#%%
fwfghads

#%%
m1 = .127
m2 = .021
x=np.linspace(0,300,300)
y1 = np.exp(-x*m1)
y2 = np.exp(-x*m2)
cc = 1/1.1
plt.plot((2-y1-y2)*cc)
mask = (df.r==16)&(df.l==df.l0)&(df.angle==0)
plt.plot(df.l[mask],(df.D/df.ctdi)[mask],'x')


#%%
#Find remaining value
t=np.array(df.final[(df.angle==0)&(df.r==16)])


#%%
masky = (df.angle==0)&(df.l==df.l0)&(df.r==16)

df['colled']=df.Dnorm/df.kcoll
df['colled'][masky]

#%%





l=np.linspace(0,1500,500)
l0 = 40
kvp = 120

t1 = Afit(cA.A,(l0,l,kvp))
plt.plot(t1)

#%%
ddf = pd.read_excel('predict.xlsx')


def get_ksize(idf):
    ksize = idf.apply(lambda x: pl.relative_dose(x['ap'], x['lr'], x['kvp'], x['filt']), axis=1)
    kangle = idf.apply(lambda x: pl.relative_dose_to_front(x['ap'], x['lr'], x['kvp'], x['filt'], x['angle']/180*np.pi), axis=1)
    return ksize,kangle



def predict(idf,ic):

    X=[idf.l0,idf.l,idf.kvp]
    idf['klength'] = fklength(ic.A,X)
    idf['kcoll'] = fkcoll(ic.A,X)
    idf['ksize'],idf.kangle = get_ksize(idf)
    idf['kall'] = ic.A[0]*idf.klength*idf.kcoll*idf.ksize*idf.kangle*idf.p
    idf['predict'] = idf.ctdi*idf.kall
    print(idf)
    return idf


ddf = predict(ddf,cA)




#%%
maskm = (df.l==df.l0)&(df.angle==0)&(df.r==16)
dfm = df[maskm]

def mfit(A,X):
    l0 = X[0]
    kvp = X[1]
    C=A[0]
    ma=A[1]
    mb = A[2]
    return (C+A[3]*kvp)*(1-np.exp(-l0*(ma+mb*kvp)))
              
Am=[1.5,0.15,-0.0007,.01]


cm=c(np.row_stack((dfm.l0,dfm.kvp)),dfm.Dnorm,sigma=dfm.Dnormerr,fitfunc = mfit,name = 'm fit',guess = Am, xlabel = 'Beam collimation (mm)')

x=np.linspace(0,40,40)
kvs = [80,100,120,140]
plt.errorbar(dfm.l,dfm.Dnorm,dfm.Dnormerr,fmt='x')
for kv in kvs:
    plt.plot(x,mfit(cm.A,[x,kv]),label=str(kv))
plt.legend()

plt.show()

df['kcoll'] = mfit(cm.A,(df.l0,df.kvp))/cm.A[0]

def Afit(A,X):
    #X:[l0,l,kvp]
    #A:[asym1,asym2,A]
    C=A[0]
    kasym = fkasym(A[1:3],X[1:])
    kcoll = fkcoll([A[3]],X[0:3])
    D = C*kcoll*kasym
    return D

AA=[1.37,0.024,-0.000035,2]