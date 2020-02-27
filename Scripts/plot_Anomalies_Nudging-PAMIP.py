"""
Script plots anomalies for multiple SC-WACCM4 experiments from PAMIP and AA
nudging from Peings et al. 2019, [GRL]

Notes
-----
    Author : Zachary Labe
    Date   : 18 February 2020
"""

### Import modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import calc_Utilities as UT
import scipy.signal as SS
import scipy.stats as sts
import read_CTLNQ as CONT
import read_ExpMonthly as NUDG
import read_ShortCoupled as COUP
import read_SIT as THICK
import read_SIC as CONC

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/Seasons/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Composites of Nudging-PAMIP - %s----' % titletime)

### Add parameters
su = [0,1,2,3,4,5]
period = 'DJF'
cps = 'yes' 
level = 'surface'
varnames = ['SLP','Z500','U700','U200','U10',
            'Z50','T2M','T700','T500','THICK']
if cps == 'none':
    runnames = [r'AA-2030',r'AA-2060',r'AA-2090',
                r'2.3--2.1',r'$\Delta$SIT',r'$\Delta$SIC']
elif cps == 'yes':
    runnames = [r'AA-2030',r'AA-2060',r'AA-2090-cps',
            r'2.3--2.1',r'$\Delta$SIT',r'$\Delta$SIC']

### Function to read in data
def readData(simu,period,vari,level,cps):
    if any([vari=='T700',vari=='T500']):
        varia = 'TEMP'
        level = 'profile'
    elif vari == 'U700':
        varia = 'U'
        level = 'profile'
    elif any([vari=='V925',vari=='V700',vari=='V1000']):
        varia = 'V'
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
        lat,lon,lev,future = NUDG.readExperi(varia,'AA','2090',level,cps)
        lat,lon,lev,historical = CONT.readControl(varia,level,cps)
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
    elif vari == 'V925':
        levq = np.where(lev == 925)[0]
        future = future[:,:,levq,:,:].squeeze()
        historical = historical[:,:,levq,:,:].squeeze()
    elif vari == 'V700':
        levq = np.where(lev == 700)[0]
        future = future[:,:,levq,:,:].squeeze()
        historical = historical[:,:,levq,:,:].squeeze()
    elif vari == 'V1000':
        levq = np.where(lev == 1000)[0]
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
        futurem = np.nanmean(future[:,-3:,:,:],axis=1)
        historicalm = np.nanmean(historical[:,-3:,:,:],axis=1)
    elif period == 'D':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,-1:,:,:],axis=1)
        historicalm = np.nanmean(historical[:,-1:,:,:],axis=1)
    elif period == 'DJF':
        print('Calculating over %s months!' % period)
        runs = [future,historical]
        var_mo = np.empty((2,historical.shape[0]-1,historical.shape[2],historical.shape[3]))
        for i in range(len(runs)):
            var_mo[i,:,:,:] = UT.calcDecJanFeb(runs[i],runs[i],lat,lon,'surface',1) 
        futurem = var_mo[0]
        historicalm = var_mo[1]
    elif period == 'JFM':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,0:3,:,:],axis=1)
        historicalm = np.nanmean(historical[:,0:3,:,:],axis=1)
    elif period == 'JF':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,0:2,:,:],axis=1)
        historicalm = np.nanmean(historical[:,0:2,:,:],axis=1)
    elif period == 'FMA':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,1:4,:,:],axis=1)
        historicalm = np.nanmean(historical[:,1:4,:,:],axis=1)
    elif period == 'FM':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,1:3,:,:],axis=1)
        historicalm = np.nanmean(historical[:,1:3,:,:],axis=1)
    elif period == 'J':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,0:1,:,:],axis=1)
        historicalm = np.nanmean(historical[:,0:1,:,:],axis=1)
    elif period == 'F':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,1:2,:,:],axis=1)
        historicalm = np.nanmean(historical[:,1:2,:,:],axis=1)
    elif period == 'M':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,2:3,:,:],axis=1)
        historicalm = np.nanmean(historical[:,2:3,:,:],axis=1)
    elif period == 'MA':
        print('Calculating over %s months!' % period)
        futurem = np.nanmean(future[:,2:4,:,:],axis=1)
        historicalm = np.nanmean(historical[:,2:4,:,:],axis=1)
    else:
        print(ValueError('Selected wrong month period!'))

    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    ### Calculate anomalies
    anom = futurem - historicalm
    
    ### Calculate ensemble mean
    anommean = np.nanmean(anom,axis=0)
    
    ### Calculate significance
    pruns = UT.calc_FDR_ttest(futurem[:,:,:],historicalm[:,:,:],0.05) #FDR
    
    ### Select climo
    climo = np.nanmean(historicalm,axis=0)
    
    return lat,lon,lev,anommean,nens,pruns,climo

############################################################################### 
############################################################################### 
############################################################################### 
### Call functions
for rr in range(len(varnames)):
    lat,lon,lev,anomAA30,nensAA30,prunsAA30,climoAA30 = readData('AA-2030',period,varnames[rr],level,cps)
    lat,lon,lev,anomAA60,nensAA60,prunsAA60,climoAA60 = readData('AA-2060',period,varnames[rr],level,cps)
    lat,lon,lev,anomAA90,nensAA90,prunsAA90,climoAA90 = readData('AA-2090',period,varnames[rr],level,cps)
    lat,lon,lev,anomcoup,nensCOUP,prunsCOUP,climoCOUP = readData('coupled',period,varnames[rr],level,cps)
    lat,lon,lev,anomthic,nensTHIC,prunsTHIC,climoTHIC = readData('SIT',period,varnames[rr],level,cps)
    lat,lon,lev,anomconc,nensCONC,prunsCONC,climoCONC = readData('SIC',period,varnames[rr],level,cps)
    
    ### Chunk data
    dataall = [anomAA30,anomAA60,anomAA90,anomcoup,anomthic,anomconc]
    nensall = [nensAA30,nensAA60,nensAA90,nensCOUP,nensTHIC,nensCONC]
    pall =    [prunsAA30,prunsAA60,prunsAA90,prunsCOUP,prunsTHIC,prunsCONC]
    climoall =[climoAA30,climoAA60,climoAA90,climoCOUP,climoTHIC,climoCONC]
         
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot anomaly composites for each experiment
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[rr] == 'T2M':
        limit = np.arange(-15,15.01,0.5)
        barlim = np.arange(-15,16,5)
        cmap = cmocean.cm.balance
        label = r'\textbf{$^{\circ}$C}'
    elif varnames[rr] == 'T700':
        limit = np.arange(-5,5.01,0.25)
        barlim = np.arange(-5,6,5)
        cmap = cmocean.cm.balance
        label = r'\textbf{$^{\circ}$C}'
    elif varnames[rr] == 'T500':
        limit = np.arange(-3,3.01,0.25)
        barlim = np.arange(-3,4,3)
        cmap = cmocean.cm.balance
        label = r'\textbf{$^{\circ}$C}'
    elif varnames[rr] == 'Z500':
        limit = np.arange(-50,50.1,1)
        barlim = np.arange(-50,51,25)
        cmap = cmocean.cm.balance
        label = r'\textbf{m}'
    elif varnames[rr] == 'Z50':
        limit = np.arange(-50,50.1,1)
        barlim = np.arange(-50,51,25)
        cmap = cmocean.cm.balance
        label = r'\textbf{m}'
    elif varnames[rr] == 'U200':
        limit = np.arange(-5,5.1,0.25)
        barlim = np.arange(-5,6,5)
        cmap = cmocean.cm.balance
        label = r'\textbf{m/s}'
    elif varnames[rr] == 'U700':
        limit = np.arange(-5,5.1,0.25)
        barlim = np.arange(-5,6,5)
        cmap = cmocean.cm.balance
        label = r'\textbf{m/s}'
    elif varnames[rr] == 'U10':
        limit = np.arange(-5,5.1,0.25)
        barlim = np.arange(-5,6,5)
        cmap = cmocean.cm.balance
        label = r'\textbf{m/s}'
    elif varnames[rr] == 'V925':
        limit = np.arange(-1,1.1,0.05)
        barlim = np.arange(-1,2,1)
        cmap = cmocean.cm.balance
        label = r'\textbf{m/s}'
    elif varnames[rr] == 'V700':
        limit = np.arange(-1,1.1,0.05)
        barlim = np.arange(-1,2,1)
        cmap = cmocean.cm.balance
        label = r'\textbf{m/s}'
    elif varnames[rr] == 'V1000':
        limit = np.arange(-1,1.1,0.05)
        barlim = np.arange(-1,2,1)
        cmap = cmocean.cm.balance
        label = r'\textbf{m/s}'
    elif varnames[rr] == 'SLP':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
        cmap = cmocean.cm.balance
        label = r'\textbf{hPa}'
    elif varnames[rr] == 'THICK':
        limit = np.arange(-40,40.1,5)
        barlim = np.arange(-40,41,20)
        cmap = cmocean.cm.balance
        label = r'\textbf{m}'
    elif varnames[rr] == 'SST':
        limit = np.arange(-1,1.01,0.05)
        barlim = np.arange(-1,2,1)
        cmap = cmocean.cm.balance
        label = r'\textbf{$^{\circ}$C}'
    elif varnames[rr] == 'RNET':
        limit = np.arange(-100,100.1,1)
        barlim = np.arange(-100,101,50)
        cmap = cmocean.cm.balance
        label = r'\textbf{W m$^{-2}$}'       
    fig = plt.figure()
    for i in range(len(runnames)):
        var = dataall[i]    
        pvar = pall[i]
        clim = climoall[i]
    
        ax1 = plt.subplot(2,3,su[i]+1)
        
        if varnames[rr] == 'SST':
            m = Basemap(projection='moll',lon_0=0,resolution='l')   
        else:
            m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                        area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        pvar, lons_cyclic = addcyclic(pvar, lon)
        pvar, lons_cyclic = shiftgrid(180., pvar, lons_cyclic, start=False)
        clim, lons_cyclic = addcyclic(clim, lon)
        clim, lons_cyclic = shiftgrid(180., clim, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
           
        circle = m.drawmapboundary(fill_color='white',color='dimgray',
                          linewidth=0.7)
        circle.set_clip_on(False)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['.....'])
        if varnames[rr] == 'Z50': # the interval is 250 m 
            cs2 = m.contour(x,y,clim,np.arange(21900,23500,250),
                            colors='k',linewidths=1.1,zorder=10)
                  
        m.drawcoastlines(color='dimgray',linewidth=0.7)
        if varnames[rr] == 'SST':
            m.fillcontinents(color='dimgray')
                
        cs.set_cmap(cmap) 
        ax1.annotate(r'\textbf{%s}' % runnames[i],xy=(0,0),xytext=(0.865,0.91),
                     textcoords='axes fraction',color='k',fontsize=11,
                     rotation=320,ha='center',va='center')
        ax1.annotate(r'\textbf{[%s]}' % nensall[i],xy=(0,0),xytext=(0.085,0.91),
                     textcoords='axes fraction',color='dimgrey',fontsize=8,
                     rotation=0,ha='center',va='center')
    
    ###########################################################################
    cbar_ax = fig.add_axes([0.293,0.1,0.4,0.03])             
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    cbar.set_label(label,fontsize=11,color='dimgrey',labelpad=1.4)  
    
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01,labelsize=8)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.16,wspace=0,hspace=0.01)
    
    if cps == 'none':
        plt.savefig(directoryfigure + 'Composites_NUDGE-PAMIP_%s_%s.png' % (period,varnames[rr]),
                    dpi=300)
    elif cps == 'yes':        
        plt.savefig(directoryfigure + 'CPS/Composites_NUDGE-PAMIP_%s_%s.png' % (period,varnames[rr]),
                    dpi=300)