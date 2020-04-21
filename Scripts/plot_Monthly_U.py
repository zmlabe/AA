"""
Script plots seasonal cycle of eddy-driven jet
Notes
-----
    Author : Zachary Labe
    Date   : 24 February 2020
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
import read_SIT_E3SM as E3SIT
import read_SIC_E3SM as E3SIC
import read_OldIceExperi as OLD
import read_LongCoupled as LC

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/SeasonalCycle/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Eddy-driven Jet %s----' % titletime)

### Add parameters
datareader = True
latpolar = 65.
variable = 'U700'
period = 'timemonth' 
level = 'surface'
runnames = [r'$\Delta$AA-2030',r'$\Delta$AA-2060',r'$\Delta$AA-2090',
            r'$\Delta$SIC-Pd',r'$\Delta$S-Coupled-Pd',r'$\Delta$SIT-Pd']
#runnamesdata = ['AA-2030','AA-2060','AA-2090','SIC','SIT','SIC']
monthstext = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m"]

### Function to read in data
def readData(simu,period,vari,level,latpolar):
    if vari == 'U700':
        varia = 'U'
        level = 'profile'
    else:
        varia = vari
    
    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    if simu == 'AA-2030':
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2030',level,'none')
        lat,lon,lev,historical = CONT.readControl(varia,level,'none')
    elif simu == 'AA-2060':
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2060',level,'none')
        lat,lon,lev,historical = CONT.readControl(varia,level,'none')
    elif simu == 'AA-2090':
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2090',level,'none')
        lat,lon,lev,historical = CONT.readControl(varia,level,'none')
    ############################################################################### 
    elif simu == 'coupled_Pd':
        lat,lon,lev,future = COUP.readCOUPs(varia,'C_Fu',level)
        lat,lon,lev,historical = COUP.readCOUPs(varia,'C_Pd',level)      
    ############################################################################### 
    elif simu == 'coupled_Pi':
        lat,lon,lev,future = COUP.readCOUPs(varia,'C_Fu',level)
        lat,lon,lev,historical = COUP.readCOUPs(varia,'C_Pi',level)  
    ###############################################################################        
    elif simu == 'SIT':
        lat,lon,lev,future = THICK.readSIT(varia,'SIT_Fu',level)
        lat,lon,lev,historical = THICK.readSIT(varia,'SIT_Pd',level)
    ############################################################################### 
    elif simu == 'SIC_Pd':
        lat,lon,lev,future = CONC.readSIC(varia,'Fu',level)
        lat,lon,lev,historical = CONC.readSIC(varia,'Pd',level)
    ############################################################################### 
    elif simu == 'SIC_Pi':
        lat,lon,lev,future = CONC.readSIC(varia,'Fu',level)
        lat,lon,lev,historical = CONC.readSIC(varia,'Pi',level)
    ############################################################################### 
    elif simu == 'E3SIT':
        lat,lon,lev,future = E3SIT.readE3SM_SIT(varia,'ESIT_Fu',level)
        lat,lon,lev,historical = E3SIT.readE3SM_SIT(varia,'ESIT_Pd',level)
    ############################################################################### 
    elif simu == 'E3SIC_Pd':
        lat,lon,lev,future = E3SIC.readE3SM_SIC(varia,'ESIC_Fu',level)
        lat,lon,lev,historical = E3SIC.readE3SM_SIC(varia,'ESIC_Pd',level)
    elif simu == 'E3SIC_Pi':
        lat,lon,lev,future = E3SIC.readE3SM_SIC(varia,'ESIC_Fu',level)
        lat,lon,lev,historical = E3SIC.readE3SM_SIC(varia,'ESIC_Pi',level)
    ############################################################################### 
    elif simu == 'OLD':
        lat,lon,lev,future = OLD.readOldIceExperi(varia,'FICT',level)
        lat,lon,lev,historical = OLD.readOldIceExperi(varia,'HIT',level)
    ############################################################################### 
    elif simu == 'LONG':
        lat,lon,lev,future = LC.readLong(varia,'Long_Fu',level)
        lat,lon,lev,historical = LC.readLong(varia,'Long_Pd',level)
    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    ### Calculate number of ensembles
    nens = np.shape(historical)[0]
    
    ### Check for 4D field
    if vari == 'T700':
        levq = np.where(lev == 700)[0]
        future = future[:,:,levq,:,:].squeeze()
        historical = historical[:,:,levq,:,:].squeeze()
    elif vari == 'T500':
        levq = np.where(lev == 500)[0]
        future = future[:,:,levq,:,:].squeeze()
        historical = historical[:,:,levq,:,:].squeeze()
    elif vari == 'U700':
        levq = np.where(lev == 700)[0]
        future = future[:,:,levq,:,:].squeeze()
        historical = historical[:,:,levq,:,:].squeeze()

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
    elif period == 'annual':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future,axis=1)
        historicalm = np.nanmean(historical,axis=1)
    elif period == 'NONE':
        print('Calculating over %s months!' % period)
        futurem = future
        historicalm = historical
    elif period == 'timemonth':
        print('Calculating over O,N,D,J,F,M months!')
        futurem = np.append(future[:,-3:,:,:],future[:,:3,:,:],axis=1)
        historicalm = np.append(historical[:,-3:,:,:],historical[:,:3,:,:],axis=1)
    else:
        print(ValueError('Selected wrong month period!'))

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Calculate polar cap
    lon2,lat2 = np.meshgrid(lon,lat)
    
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
lat,lon,lev,anomAA30,nensAA30,prunsAA30,climoAA30 = readData('AA-2030',period,variable,level,latpolar)
lat,lon,lev,anomAA60,nensAA60,prunsAA60,climoAA60 = readData('AA-2060',period,variable,level,latpolar)
lat,lon,lev,anomAA90,nensAA90,prunsAA90,climoAA90 = readData('AA-2090',period,variable,level,latpolar)
lat,lon,lev,anomcoup,nensCOUP,prunsCOUP,climoCOUP = readData('SIC_Pd',period,variable,level,latpolar)
lat,lon,lev,anomthic,nensTHIC,prunsTHIC,climoTHIC = readData('coupled_Pd',period,variable,level,latpolar)
lat,lon,lev,anomconc,nensCONC,prunsCONC,climoCONC = readData('SIT',period,variable,level,latpolar)

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
if variable == 'U700':
    limit = np.arange(-1.5,1.51,0.1)
    limitc = np.arange(3,71,3)
    barlim = np.arange(-1.5,1.6,1.5)
    cmap = cmocean.cm.balance
    label = r'\textbf{U700 [m/s]}'
    zscale = np.arange(-90,91,15)
    time = np.arange(0,6,1)
    latq,timeq = np.meshgrid(lat,time)
elif variable == 'U200':
    limit = np.arange(-3,3.1,0.1)
    limitc = np.arange(0,71,10)
    barlim = np.arange(-3,4,3)
    cmap = cmocean.cm.balance
    label = r'\textbf{U200 [m/s]}'
    zscale = np.arange(-90,91,15)
    time = np.arange(0,6,1)
    latq,timeq = np.meshgrid(lat,time)
elif variable == 'U10':
    limit = np.arange(-5,5.1,0.1)
    limitc = np.arange(-70,71,5)
    barlim = np.arange(-5,6,5)
    cmap = cmocean.cm.balance
    label = r'\textbf{U10 [m/s]}'
    zscale = np.arange(-90,91,15)
    time = np.arange(0,6,1)
    latq,timeq = np.meshgrid(lat,time)
        
fig = plt.figure()
for i in range(len(runnames)):
    
    var = dataall[i]
    pvar = pall[i]
    clim = climoall[i]
    en = nensall[i]
    
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
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',color='k',fontsize=7)
    elif i == 3:
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')   
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
        plt.gca().axes.get_xaxis().set_visible(True)
        plt.gca().axes.get_yaxis().set_visible(True)
        plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',color='k',fontsize=7)
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
    
    ### Plot contours
    cs = plt.contourf(timeq,latq,var,limit,extend='both')
    cs1 = plt.contour(timeq,latq,clim,limitc,colors='dimgrey',
                      linewidths=1)
    cs2 = plt.contourf(timeq,latq,pvar,colors='None',
                   hatches=['//////'],linewidths=0.4)
    cs.set_cmap(cmap)
    
    plt.xticks(np.arange(0,6,1),monthstext,fontsize=4)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    
    plt.xlim([0,5])
    plt.ylim([0,90])
    plt.minorticks_off()
           
    ax1.annotate(r'\textbf{%s}' % runnames[i],xy=(0,90),xytext=(0.98,0.93),
         textcoords='axes fraction',color='k',fontsize=8,
         rotation=0,ha='right',va='center')
    ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,90),xytext=(0.02,0.93),
         textcoords='axes fraction',color='k',fontsize=8,
         rotation=0,ha='left',va='center')
    ax1.annotate(r'\textbf{[%s]}' % en,xy=(0,90),xytext=(0.02,0.07),
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
plt.savefig(directoryfigure + 'MonthlyModels_%s_Oct-Apr.png' % variable,dpi=300)
print('Completed: Script done!')