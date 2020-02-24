"""
Script plots relationship between vertical warming in modeled data sets

Notes
-----
    Author : Zachary Labe
    Date   : 21 February 2020
"""

### Import modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import calc_Utilities as UT
import scipy.stats as sts
import read_CTLNQ as CONT
import read_ExpMonthly as NUDG
import read_ShortCoupled as COUP
import read_SIT as THICK
import read_SIC as CONC

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/Vertical_Model/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Vertical Warming- %s----' % titletime)

### Add parameters
datareader = True
latpolar = 65.
variable = 'TEMP'
period = 'FM' 
level = 'profile'
runnames = [r'AA-2030',r'AA-2060',r'AA-2090',
            r'2.3--2.1',r'$\Delta$SIT',r'$\Delta$SIC']
runnamesdata = ['AA-2030','AA-2060','AA-2090','coupled','SIT','SIC']

### Function to read in data
def readData(simu,period,varia,level):
    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    if simu == 'AA-2030':
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2030',level)
        lat,lon,lev,historical = CONT.readControl(varia,level)
    elif simu == 'AA-2060':
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2060',level)
        lat,lon,lev,historical = CONT.readControl(varia,level)
    elif simu == 'AA-2090':
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2090',level)
        lat,lon,lev,historical = CONT.readControl(varia,level)
    ############################################################################### 
    elif simu == 'coupled':
        lat,lon,lev,future = COUP.readCOUPs(varia,'C_Fu',level)
        lat,lon,lev,historical = COUP.readCOUPs(varia,'C_Pd',level)        
    ###############################################################################        
    elif simu == 'SIT':
        lat,lon,lev,future = THICK.readSIT(varia,'SIT_Fu',level)
        lat,lon,lev,historical = THICK.readSIT(varia,'SIT_Pd',level)
    ############################################################################### 
    elif simu == 'SIC':
        lat,lon,lev,future = CONC.readSIC(varia,'Fu',level)
        lat,lon,lev,historical = CONC.readSIC(varia,'Pd',level)
    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    ### Calculate number of ensembles
    nens = np.shape(historical)[0]

    ### Check for missing data [ensembles,months,lat,lon]
    future[np.where(future <= -1e10)] = np.nan
    historical[np.where(historical <= -1e10)] = np.nan
    
    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    ### Calculate over period
    if period == 'OND':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,-3:],axis=1)
        historicalm = np.nanmean(historical[:,-3:],axis=1)
    elif period == 'D':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,-1:],axis=1)
        historicalm = np.nanmean(historical[:,-1:],axis=1)
    elif period == 'DJF':
        print('Calculating over %s months!' % period)
        runs = [future,historical]
        var_mo = np.empty((2,historical.shape[0]-1,historical.shape[2],historical.shape[3],historical.shape[4]))
        for i in range(len(runs)):
            var_mo[i,:,:,:,:] = UT.calcDecJanFeb(runs[i],runs[i],lat,lon,level,17) 
        futurem = var_mo[0]
        historicalm = var_mo[1]
    elif period == 'JFM':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,0:3],axis=1)
        historicalm = np.nanmean(historical[:,0:3],axis=1)
    elif period == 'JF':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,0:2],axis=1)
        historicalm = np.nanmean(historical[:,0:2],axis=1)
    elif period == 'FMA':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,1:4],axis=1)
        historicalm = np.nanmean(historical[:,1:4],axis=1)
    elif period == 'FM':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,1:3],axis=1)
        historicalm = np.nanmean(historical[:,1:3],axis=1)
    elif period == 'J':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,0:1],axis=1)
        historicalm = np.nanmean(historical[:,0:1],axis=1)
    elif period == 'F':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,1:2],axis=1)
        historicalm = np.nanmean(historical[:,1:2],axis=1)
    elif period == 'M':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,2:3],axis=1)
        historicalm = np.nanmean(historical[:,2:3],axis=1)
    elif period == 'MA':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,2:4],axis=1)
        historicalm = np.nanmean(historical[:,2:4],axis=1)
    else:
        print(ValueError('Selected wrong month period!'))

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Calculate zonal means
    futuremz = np.nanmean(futurem,axis=3)
    historicalmz = np.nanmean(historicalm,axis=3)
    
    ### Calculate anomalies [ens,level,lat]
    anom = futuremz - historicalmz

    ### Calculate ensemble mean
    anommean = np.nanmean(anom,axis=0)
    
    ### Calculate significance
    pruns = UT.calc_FDR_ttest(futuremz,historicalmz,0.05) #FDR
    
    ### Select climo
    climo = np.nanmean(historicalmz,axis=0)
    
    return lat,lon,lev,anommean,nens,pruns,climo

### Call data
lat,lon,lev,anomAA30,nensAA30,prunsAA30,climoAA30 = readData('AA-2030',period,variable,level)
lat,lon,lev,anomAA60,nensAA60,prunsAA60,climoAA60 = readData('AA-2060',period,variable,level)
lat,lon,lev,anomAA90,nensAA90,prunsAA90,climoAA90 = readData('AA-2090',period,variable,level)
lat,lon,lev,anomcoup,nensCOUP,prunsCOUP,climoCOUP = readData('coupled',period,variable,level)
lat,lon,lev,anomthic,nensTHIC,prunsTHIC,climoTHIC = readData('SIT',period,variable,level)
lat,lon,lev,anomconc,nensCONC,prunsCONC,climoCONC = readData('SIC',period,variable,level)

### Chunk data
dataall = [anomAA30,anomAA60,anomAA90,anomcoup,anomthic,anomconc]
nensall = [nensAA30,nensAA60,nensAA90,nensCOUP,nensTHIC,nensCONC]
pall =    [prunsAA30,prunsAA60,prunsAA90,prunsCOUP,prunsTHIC,prunsCONC]
climoall =[climoAA30,climoAA60,climoAA90,climoCOUP,climoTHIC,climoCONC]

###########################################################################
###########################################################################
###########################################################################
##### Plot profiles
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
        
### Set limits for contours and colorbars
if variable == 'TEMP':
    limit = np.arange(-3,3.01,0.25)
    barlim = np.arange(-3,4,1)
    cmap = cmocean.cm.balance
    label = r'\textbf{$^{\circ}$C}'
    zscale = np.array([1000,925,850,700,500,300,200])
    latq,levq = np.meshgrid(lat,lev)
        
fig = plt.figure()
for i in range(len(runnames)):
    
    varnomask = dataall[i]
    pvar = pall[i]
    clim = climoall[i]
    en = nensall[i]
    
    ### Mask significant
    pvar[np.isnan(pvar)] = 0.
    var = varnomask * pvar
    var[var == 0.] = np.nan
    
    ### Create plot
    ax1 = plt.subplot(2,3,i+1)
    ax1.spines['top'].set_color('dimgrey')
    ax1.spines['right'].set_color('dimgrey')
    ax1.spines['bottom'].set_color('dimgrey')
    ax1.spines['left'].set_color('dimgrey')
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['right'].set_linewidth(2)
    ax1.spines['top'].set_linewidth(2)
    if i == 0:
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
        plt.gca().axes.get_yaxis().set_visible(True)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.ylabel(r'\textbf{Pressure [hPa]}',color='k',fontsize=7)
    elif i == 3:
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')   
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
        plt.gca().axes.get_xaxis().set_visible(True)
        plt.gca().axes.get_yaxis().set_visible(True)
        plt.ylabel(r'\textbf{Pressure [hPa]}',color='k',fontsize=7)
    elif i == 4 or i == 5:
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')   
        plt.gca().axes.get_xaxis().set_visible(True)
        plt.gca().axes.get_yaxis().set_visible(False)
    else:
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
            width=0,color='w')
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.gca().axes.get_xaxis().set_visible(False)

    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    
    ### Add levels
    plt.axhline(lev[1],linewidth=0.6,linestyle='--',dashes=(1,0.5),
                color='k')
    plt.axhline(lev[2],linewidth=0.6,linestyle='--',dashes=(1,0.5),
                color='k')
    plt.axhline(lev[3],linewidth=0.6,linestyle='--',dashes=(1,0.5),
                color='k')
    plt.axhline(lev[5],linewidth=0.6,linestyle='--',dashes=(1,0.5),
                color='k')
    plt.axhline(lev[7],linewidth=0.6,linestyle='--',dashes=(1,0.5),
                color='k')
    
    ### Plot contours
    cs = plt.contourf(latq,levq,var,limit,extend='both')
#    cs1 = plt.contour(latq,levq,clim,np.arange(-90,95,5),
#                      linewidths=0.3,colors='k',extend='both')
#    cs2 = plt.contourf(latq,levq,pvar,colors='None',
#                   hatches=['//////'],linewidths=0.4)
    cs.set_cmap(cmap) 
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=6)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    
    plt.xlim([45,90])
    if variable == 'TEMP' or variable == 'GEOP':
        plt.ylim([1000,200])
    else:
        plt.ylim([1000,10])
    plt.minorticks_off()
           
    ax1.annotate(r'\textbf{%s}' % runnames[i],xy=(80,200),xytext=(0.98,0.93),
         textcoords='axes fraction',color='k',fontsize=8,
         rotation=0,ha='right',va='center')
    ax1.annotate(r'\textbf{[%s]}' % en,xy=(80,200),xytext=(0.02,0.93),
         textcoords='axes fraction',color='dimgrey',fontsize=8,
         rotation=0,ha='left',va='center')

###########################################################################
plt.tight_layout()
cbar_ax = fig.add_axes([0.33,0.08,0.4,0.03])                    
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(label,fontsize=11,color='dimgrey',labelpad=1.4)  

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.001,labelsize=7)
cbar.outline.set_edgecolor('dimgrey')
    
plt.subplots_adjust(bottom=0.17,hspace=0.08,wspace=0.08)    
plt.savefig(directoryfigure + 'VerticalModels_TEMP_%s.png' % period,
            dpi=300)
print('Completed: Script done!')