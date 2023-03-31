import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc
from shapely.geometry import LineString
from shapely.geometry import Point
import shapely.affinity
from scipy.signal import savgol_filter

import os
script_dir = os.path.dirname(__file__)


R=55
#%%
class shape():
    #elipse definition
    def __init__(self,r_AP,r_LR):
        self.AP = r_AP
        self.LR = r_LR

    
    def xy(self,theta):
        return self.LR*np.sin(theta),self.AP*np.cos(theta)
        
    def r(self,theta):
        x,y = self.xy(theta)
        return np.sqrt(x**2+y**2)
    


#%%
def get_geo(AP,LR,angle=0,CT_rad = R):
    
    ellipse = shape(AP,LR)
        
    P = np.array(ellipse.xy(angle))[:,np.newaxis]
    
    #Set number of angles to be considered. 
    n_angles = 360
    #Construct an input variable for angle, in radians
    theta = np.linspace(0,360,n_angles+1)
    thetar = theta/180*np.pi
    
    #Vector from isocentre to source at each angle
    Q = np.array([R*np.sin(thetar),CT_rad*np.cos(thetar)])
    
    #Vector from initial point to source at each angle
    D=Q-P
    #Distance from initial point to source at each angle
    Dlength = np.sqrt(D[0,:]**2 + D[1,:]**2)
    
    #Dot product for calculating angle between point and beam midline
    DdotQ = np.einsum('ij,ij->j',Q,D)
    phi = np.arccos(DdotQ/Dlength/CT_rad)
    
    
    
    #Figure out the path length through the phantom
    #Use a line consisting of n_points
    n_points = 1001
    
    #Line runs along D, vector from initial point to source at each angle
    #Constructing the line vector:
    #Make a linspace
    test = np.linspace(0,1,n_points)
    test = test[:,np.newaxis]
    
    #This is a fancey eigensum. Create a 3d array:
    #axis 1: distance along line, out of n_poitns
    #axis 2: gives x,y coordinates for vector
    #axis 3: for each angle 
    #n.b. this is fucking sorcery, you're welcome future me
    testvectors = np.einsum('ij,jk->ijk',test,D)
    
    #Add P to get absolute vector to each point along the line
    testlines = testvectors+P
    
    #Set any 0s to a small number, to eliminate divide by 0 issue
    testlines[testlines==0] = 0.000001
    
    #Pythagorean sum of vectors to get from shape centre to line point
    testlines_rad = np.sqrt(testlines[:,0,:]**2+testlines[:,1,:]**2)
    #cosine to get angle between isocentre and line point for each line point
    #Fuck you past me for not commenting this sorcery
    testlines_angle = np.arctan(testlines[:,0,:]/testlines[:,1,:])
    #Get radius of shape for each and every point
    shape_rad = ellipse.r(testlines_angle)

    #make binary mask for each point on each line that is within shape radius    
    binary = testlines_rad<=shape_rad
    #Total number of points on each line that are within shape radius
    dlength = np.sum(binary,axis=0)
    #Multiply number of segments within radius by length of each segment to get dlength
    segmentlength = Dlength/n_points
    dlength = dlength*segmentlength

    return (theta, Dlength,dlength,phi),ellipse

#%%
#Load depth dose data

#The depth dose class. Initialise with 
class DD():
    def __init__(self,fn):
        df = pd.read_excel(fn)
        df.Y = savgol_filter(df.Y,101,1)
        maxindex=df[df.Y == df.Y.max()].index[0]
    
        D = df.Y
        #[maxindex:
        d = np.arange(0,len(D))/50-maxindex/50
        D = D * (d+30)**2
        D=D/D.max()
        maxindex=D[D == D.max()].index[0]
        D=D[d<15][maxindex:]
        d=d[d<15][maxindex:]
        self.d=d
        self.D = D
        self.fit_cutoff = 0.7
        self.dose_cutoff = 0.9
        self.fit_exp()
    
    def fit_exp(self):
        cutoff_mask = self.d > np.max(self.d)*self.fit_cutoff
        plt.show()
        self.m,self.A = np.polyfit(self.d[cutoff_mask], np.log(self.D[cutoff_mask]), 1)
        self.A = np.exp(self.A)
        
        
    def get_D(self,d):
        d = np.array(d)
        #Make mask showing which lengths are below about 10
        mask = d<=np.max(self.d)*self.dose_cutoff
        #Make an empty array
        dose = np.zeros(d.shape)
        #For lengths below ~10, use interpolation
        dose[mask] = np.interp(d[mask],self.d,self.D)
        #Otherwise use fitted exponential. Allows depths greater than
        dose[~mask] = self.A*np.exp(self.m*d[~mask])
        return dose
    
    def get_d(self,D):
        D=np.array(D)
        Dcutoff = self.get_D(np.max(self.d)*self.dose_cutoff)
        mask = D>Dcutoff
        d = np.zeros(D.shape)
        d[mask]=np.interp(D[mask],self.D[::-1],self.d[::-1])
        d[~mask] = np.log(D[~mask]/self.A)/self.m
#        plt.plot(d)
        return d
            
        
class bowtie():
    def __init__(self,fn):
        BTdata = np.loadtxt(fn,skiprows=0,delimiter=',')
        self.x = BTdata[:,0]
        self.I = BTdata[:,1]
        self.I = savgol_filter(self.I, 151, 1)
#        self.I = self.I*.2/self.I.min()
        mask = ~np.isnan(self.I)
        self.I = self.I[mask]
        self.x = self.x[mask]
        self.phi = np.arcsin(self.x/R)
        #self.I = self.I/np.max(self.I)
        #self.d = DD.get_d(self.I)

    
    def get_I(self,phi):
        return np.interp(phi,self.phi,self.I)

    

class table():
    def __init__(self,ishape,fn):
        #read file
        #set attenuation
        #set table size/shape?
        self.shape = ishape
        self.att = 0.9
        self.width = 60 #cm
        
    def intersection(self,P_angle):
        P1 = self.shape.xy(P_angle)
        TE = np.array([self.width/2,self.shape.xy(np.pi)[1]]) #table edge point
        l = LineString([P1,P1+(TE-P1)*10]) #Extrapolate line by factor of 10 in direction of table edge
        c = Point(0,0).buffer(R).boundary #type(circle)=polygon
        mp = c.intersection(l)
        if mp.is_empty:
            return np.array(180)
        elif mp.geom_type == 'Point':
            return np.arctan(mp.x/mp.y)
        else:
            raise ValueError('something unexpected: ' + mp.geom_type)
        
    
#table(shape(16,16),'test').intersection(0)*180/np.pi

            
            


        
def dose_series(iDD,iBT,ishape,theta,Dlength,dlength,phi):
    #Load bowtie data
    bt = iBT.get_I(phi)
#    plt.plot(bt)
#    plt.show()

    #Find angle subtended by table
    itable = table(ishape,'table.csv')
    #Add effective path length through water equivalent to the table thickness
    table_angle = itable.intersection(0)
    tablemask = (theta*np.pi/180>np.pi+table_angle) & (theta*np.pi/180<np.pi-table_angle)
    table_attenuate = np.ones(phi.shape)
    table_attenuate[tablemask] = itable.att

    dose = 1/Dlength**2*iDD.get_D(dlength)*bt*table_attenuate
    return dose
    

#Load depth dose and bowtie data






    

class device_data():
    def __init__(self,dd_list,bt_list):
        self.set_dd_list(dd_list)
        self.set_bt_list(bt_list)
            
    def set_dd_list(self,dd_list):
        self.dd={}
        for dd in dd_list:
            self.set_dd(*dd)
            
    def set_bt_list(self,bt_list):
        self.bt = {}
        for bt in bt_list:
            self.set_bt(*bt)
    
    def set_dd(self,kv,fn):
        self.dd[kv] = DD(fn)
    
    def set_bt(self,filter_type,kv,fn):
        if filter_type not in self.bt:
            self.bt[filter_type] = {}
        self.bt[filter_type][kv]=bowtie(fn)
    
Devices = device_data([],[])

default_dds = [[80,'dat/dd/100.xlsx'],
                [100,'dat/dd/100.xlsx'],
                [120,'dat/dd/120.xlsx'],
                [140,'dat/dd/140.xlsx']]

def default_bts_maker():
    defaulty = []
    kvs = [80,100,120,140]
    types = ['body','head']
    bt_dir = 'dat/bowtie/'
    for kv in kvs:
        for filt in types:
            defaulty.append([filt,kv,bt_dir+filt+'.csv'])
    return defaulty
default_bts = default_bts_maker()

Devices.set_dd_list(default_dds)
Devices.set_bt_list(default_bts)
  
    
def plot_bowties(Devices=Devices):
#    fig, axes = plt.subplots(nrows=1, ncols=2,sharey=True, figsize=(9, 4))# , sharex=True, sharey=True)
    fig,ax = plt.subplots()
    for i,btname in enumerate(Devices.bt):
        ax.plot(Devices.bt[btname][120].phi*180/np.pi,Devices.bt[btname][120].I,label = btname.capitalize()+' filter')
    ax.set_xlabel(r'Angular offset from central axis ($\phi^{\circ}$)')
    ax.set_ylabel('Relative intensity')
    ax.legend(loc=1)
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    fig.tight_layout()
    fig.savefig('out/bowties.eps',format='eps',dpi=600)
    plt.show()
#plot_bowties()

def plot_dds(Devices=Devices):
    fig,ax = plt.subplots()
    for i,ddname in enumerate(Devices.dd):
        ax.plot(Devices.dd[ddname].d,Devices.dd[ddname].get_D(Devices.dd[ddname].d),label = str(ddname)+' kVp')
#        ax.plot(Devices.dd[ddname].d,Devices.dd[ddname].D,label = str(ddname)+' kVp')
    ax.set_xlabel(r'Depth in water (cm))')
    ax.set_ylabel('Relative dose')
    ax.legend(loc=1)
#    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_ybound(lower = 0)
    fig.tight_layout()
    fig.savefig('out/dds.svg',format='svg',dpi=600)
    plt.show()
    

    
    
#plot_dds()


#%%
#Get total dose for each ellipse, for circles
def total_dose(AP,LR,ddname,filter_type,angle = 0):
    try:
        iDD = Devices.dd[ddname]
    except KeyError:
        print('Could not load depth dose with name' + str(ddname))
        print('Available types are:')
        print(list(dd.keys()))
        try:
            iDD = Devices.dd['120']
            print('Defaulting to 120 kVp')
        except:
            iDD = Devices.dd[list(Devices.dd.keys())[-1]]
            print('Defaulting to '+ list(Devices.dd.keys())[-1] +'kVp')
    try:
        iBT = Devices.bt[filter_type][ddname]
    except KeyError:
        print('Could not load bowtie filter with name' + filter_type)
        print('Available types are:')
        print(list(Devices.bt.keys()))
        print('Defaulting to '+list(Devices.bt.keys())[-1])
        iBT = Devices.bt[list(Devices.bt.keys())[-1]]
    data,ishape = get_geo(AP,LR,angle)
    dose = dose_series(iDD,iBT,ishape,*data)
    return np.sum(dose)

def total_dose_series(APs,LRs,ddname,filter_type):
    totdose = [total_dose(APs[i],LRs[i],ddname,filter_type) for i in np.arange(len(APs))]
    return totdose


def relative_dose(AP,LR,ddname,filter_type):
    if (AP ==16) & (LR==16):
        return 1
    dose = total_dose(AP,LR,ddname,filter_type)
    dose16 = total_dose(16,16,ddname,'body')
    return dose/dose16

def relative_dose_to_front(AP,LR,ddname,filter_type,angle):
    return total_dose(AP,LR,ddname,filter_type,angle)/total_dose(AP,LR,ddname,filter_type,0)






#%%



def plot_geo(name,theta,D,d,phi,ax=None):
    ax.plot(theta,D,label = 'D')
    ax.plot(theta,d,label = 'd')
    ax.plot(theta,phi*180/np.pi,label = r'$\phi$' )
#    ax.set_ylabel(r'Angle ($^{\circ}$),distance (cm)')
#    ax.set_xlabel(r'$\Theta$ ($^{\circ}$)')
    ax.legend(loc=2)
    #plt.savefig('out/'+name+'.eps', format='eps', dpi=600)
    #plt.show()
    
def plot_dose(name,theta,dose,ax):
    ax.plot(theta,dose,label = 'Dose')
#    ax.set_ylabel(r'Dose (relative)')
#    ax.set_xlabel(r'$\Theta$ ($^{\circ}$)')
    #plt.legend(loc=2)
    
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    #plt.savefig('out/'+name+'_dose.eps', format='eps', dpi=600)
    #plt.show()
    





#%%


#Run for a few set shapes 
def make_geo_plots():
    geos = [(16,16),(8,8),(8,16),(16,8)]
    plot_types = ['Geometric parameters','Dose contribution']
    shape_types = ['16 cm cylinder','8 cm cylinder','8 by 16 cm ellipse','16 by 8 cm ellipse']
    
    rows = len(geos)
    #rc('text', usetex=True)
    fig, axes = plt.subplots(nrows=rows, ncols=2,sharex=True, figsize=(9, 9))# , sharex=True, sharey=True)

    for i,rows in enumerate(axes):
        AP = geos[i][0]
        LR = geos[i][1]
        name = str(AP)+'_'+str(LR)
        data,ishape = get_geo(AP,LR)
        dose = dose_series(Devices.dd[120],Devices.bt['body'][120],ishape,*data)
        
        plot_geo(name,*data,axes[i][0])
        plot_dose(name,data[0],dose,axes[i][1])
        

        for j,ax in enumerate(rows):
            if i==0:
                ax.set_title(r'\Large{{{}}}'.format(plot_types[j]))
            if i == len(axes) - 1:
                ax.set_xlabel(r'Source angle ($^{\circ}$)')
            if j == 0:
                ax.set_ylabel(r'\Large{{{}}}'.format(shape_types[i])+'\n'+r'Angle ($^{\circ}$),distance (cm)')
            else:
                ax.set_ylabel('Relative dose')

    plt.tight_layout()
    plt.savefig('out/geos_doses.svg',format='svg',dpi=600)
    #plt.show()
    #rc('text', usetex=False)
    
        
    
    
make_geo_plots()

#%%
def make_kvp_plots():
    geoAP = np.linspace(1,50,50)
    geoLR = np.linspace(1,50,50)

    for ddname in Devices.dd.keys():
        doses = total_dose_series(geoAP,geoLR,ddname,'body')/total_dose(16,16,ddname,'body')
        doses2 = total_dose_series(geoAP,geoAP*1.45,ddname,'body')/total_dose(16,16,ddname,'body')
        fig, ax = plt.subplots()
        ax.plot(geoAP*2,doses,label='Cylindrical phantom')
        ax.plot(geoAP*2,doses2,label = 'Ellipsoidal phantom')
#        ax.yaxis.set_major_formatter(plt.NullFormatter())
        ax.set_ylabel(r'Size correction factor ($k_{size}$)')
        ax.set_xlabel('Phantom AP diameter (cm)')
        ax.set_ylim(ymin=0)
        ax.axis([0,100,0.7,2])
        ax.legend(loc=2)
    
        #Show background objects
        add_patient_size(ax)
        fig.tight_layout()
        fig.savefig('out/totdose_'+str(ddname)+'.svg',format='svg',dpi=600)
        fig.show()
        

def add_patient_size(ax):
    #Head section
    ax.axvspan(15, 20, alpha=0.2, color='green')
    ax.annotate('Head',
                xy = (18,0.8),
                xytext = (3,0.8),
                arrowprops=dict(arrowstyle="->"),
                horizontalalignment='left',
                verticalalignment='center',
                                )
    #Body section
    ax.axvspan(24, 50, alpha=0.2, color='green')
    ax.annotate('Typical torso AP diameters',
                xy = (35,0.8),
                xytext = (60,0.8),
                arrowprops=dict(arrowstyle="->"),
                horizontalalignment='left',
                verticalalignment='center',
                                )

     

def make_shape_plots():
    geoAP = np.linspace(1,50,50)
    geoLR = np.linspace(1,50,50)
    shapes = [(1,1,'Cylinder'),(1,1.45,'Ellipsoid')]
    dose_list = {}
    for shape in shapes:
        dose_list[shape[2]] = {}
        for ddname in Devices.dd.keys():
            doses = total_dose_series(geoAP*shape[0],geoLR*shape[1],ddname,'body')
            dose_list[shape[2]][ddname] =doses
            plt.plot(geoAP*2,doses,label=str(ddname) + ' kVp')

        ax = plt.axes()
        ax.yaxis.set_major_formatter(plt.NullFormatter())
        plt.ylabel('Relative dose')
        plt.xlabel('Phantom AP diameter (cm)')
        plt.ylim(ymin=0)
        plt.legend(loc=2)
        
        plt.tight_layout()
        plt.savefig('out/totdose_'+shape[2]+'.svg',format='svg',dpi=600)
        plt.show()
        
    return dose_list

def make_shape_ratio_plots(data = False):
    if not data:
        data = make_shape_plots()
    a = np.array(pd.DataFrame.from_dict(data))
    a=np.array(a.tolist())
    
    print(((a[:,0,:]/a[:,1,:]).T).shape)
    plt.plot(np.arange(0,100,2),(a[:,0,:]/a[:,1,:]).T)
    plt.ylabel('Surface intensity ratio')
    plt.xlabel('Major axis phantom diameter')
    plt.legend(['80 kVp','100 kVp', '120 kVp'],loc = 1)
    plt.tight_layout()
    plt.savefig('out/shaperatio.svg',format='svg',dpi=600)
    
    return data

#make_kvp_plots()
make_shape_plots()
#%%


#t,y = get_geo(16,16)
#fig,ax = plt.subplots()
#plot_geo('test',*t,ax=ax)
#
#ax.set_xlabel(r'Source angle ($^{\circ}$)')
#ax.set_ylabel(r'Angle ($^{\circ}$),distance (cm)')
#fig.tight_layout()
#fig.savefig('out/16_16.eps',format='eps',dpi=600)
##%%
#
#dose = dose_series(Devices.dd[120],Devices.bt['body'][120],y,*t)
#
#
#fig,ax = plt.subplots()
#plot_dose('test',t[0],dose,ax)
#ax.set_xlabel(r'Source angle ($^{\circ}$)')
#ax.set_ylabel(r'Relative dose')
#fig.tight_layout()
#fig.savefig('out/16_16_dose.eps',format='eps',dpi=600)
#
#
#


#%%
#aps = np.linspace(1,45,45)
#lrs = np.linspace(1,45,45)
#
#apg,lrg = np.meshgrid(aps,lrs)
#
#def testf(ap,lr,kvp,filter):
#    return relative_dose(ap,lr,kvp,filter)
#testv = np.vectorize(testf)
#
#
#
#kvps = [80,100,120,140]
#filters = ['head','body']
#
#out = {}
#for kv in kvps:
#    out[kv]={}
#
#for kv in kvps:
#    out[kv]['body'] = testv(apg,lrg,kv,'body')
#for kv in kvps:
#    out[kv]['head'] = testv(apg,lrg,kv,'head')
#    
##%%
#    
#xl = pd.ExcelWriter('pldata.xlsx')
#
#for kv in kvps:
#    for fil in filters:
#        t = pd.DataFrame(out[kv][fil])
#        t.to_excel(xl,
#                            sheet_name = str(kv)+str(fil),
#                            header = False,
#                            index = False
#                            )
#    
#    
#    
#    
#    
#xl.save()
#    
#    
###%%
##
##from mpl_toolkits.mplot3d import Axes3D
###%%
##
##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
##
##ax.plot_surface(apg, lrg, t)
##
##ax.set_xlabel('AP Label')
##ax.set_ylabel('LR Label')
##ax.set_zlabel('ksize Label')
##
##plt.show()
#
##
##%%
##
#crap = [out[80][i,i] for i in np.arange(30)]
##%%
#
#i=16
#plt.plot(apg[i,:],t[i,:])