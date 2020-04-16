"""
Script calculates the mean vertical temperature response over the polar
cap for December-January-February for each model experiment. See text file.

Notes
-----
    Author : Zachary Labe
    Date   : 16 April 2020
"""

### Import modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
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
directoryfigure = '/home/zlabe/Desktop/AA/Emergent/Paper_v1/'
directorydata = '/home/zlabe/Documents/Research/AA/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Scatter of Warming-High (ALL)- %s----' % titletime)

### Add parameters
datareader = False
latpolar = 65.
variable = 'TEMP'
period = 'DJF' 
level = 'surface'
runnamesdata = ['AA-2030','AA-2060','AA-2090',
                'coupled_Pd','coupled_Pi','LONG','SIT',
                'SIC_Pd','SIC_Pi','OLD']
runnamesdata_E3SM = ['E3SIC_Pi','E3SIC_Pd']
runnamesdata_AMIP = ['AMIP-HL','AMIP']

def PolarCapVert(simu,varia,level,latpolar,period,levelVert):
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
        lat,lon,lev,historical = E3SIT.readE3SM_SIT(varia,'ESIT_Pd_B',level)
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
        var_mo = np.empty((2,historical.shape[0]-1,historical.shape[2],
                           historical.shape[3],historical.shape[4]))
        for i in range(len(runs)):
            var_mo[i,:,:,:,:] = UT.calcDecJanFeb(runs[i],runs[i],lat,lon,level,17) 
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
    elif period == 'NONE':
        futurem = future
        historicalm = historical
    else:
        print(ValueError('Selected wrong month period!'))

    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    ### Calculate anomalies
    anom = futurem - historicalm
    
    ### Meshgrid for lat,lon
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate SHI
    latq = np.where((lat >= latpolar))[0]
    anomp = anom[:,:,latq,:]
    lat2p = lat2[latq,:]
    polarave = UT.calc_weightedAve(anomp,lat2p)
    
    ### Calculate ensemble mean
    polaraveMean = np.nanmean(polarave,axis=0)
    
    ############################################################################### 
    ############################################################################### 
    ############################################################################### 
    ### Calculate vertical levels and save file
    levqq = np.where((lev >= levelVert))[0]
    levvv = lev[levqq]
    polaraveMeanvvv = polaraveMean[levqq]
    
    ### Save file
    np.savetxt(directorydata + '%s_1000-%s_%s.txt' % (simu,levelVert,varia),
               polaraveMeanvvv,delimiter=',',fmt='%.3f')
    np.savetxt(directorydata + 'Levels_1000-%s_%s.txt' % (levelVert,varia),
               levvv,delimiter=',',fmt='%.1f')
    
    print('\n========Calculated Polar Cap Average========\n')
    return polaraveMeanvvv,levvv

###############################################################################
###############################################################################
###############################################################################
### Test functions (do not use!)
    
#aveAA30,lev = PolarCapVert('AA-2030','TEMP','profile',65,'DJF',500)
#aveAA60,lev = PolarCapVert('AA-2060','TEMP','profile',65,'DJF',500)
#aveAA90,lev = PolarCapVert('AA-2090','TEMP','profile',65,'DJF',500)
################################################################################
#ave_coupPd,lev = PolarCapVert('coupled_Pd','TEMP','profile',65,'DJF',500)
#ave_coupPi,lev = PolarCapVert('coupled_Pi','TEMP','profile',65,'DJF',500)
################################################################################
#ave_SIT,lev = PolarCapVert('SIT','TEMP','profile',65,'DJF',500)
#ave_SICPd,lev = PolarCapVert('SIC_Pd','TEMP','profile',65,'DJF',500)
#ave_SICPi,lev = PolarCapVert('SIC_Pi','TEMP','profile',65,'DJF',500)
###############################################################################
#ave_NET,lev = PolarCapVert('OLD','TEMP','profile',65,'DJF',500)
###############################################################################
#ave_ESICPd,lev = PolarCapVert('E3SIC_Pd','TEMP','profile',65,'DJF',500)
#ave_ESICPi,lev = PolarCapVert('E3SIC_Pi','TEMP','profile',65,'DJF',500)