def PolarCap(simu,vari,level,latpolar,period):
    """
    Script calculates average over the polar cap (set latitude)
    """
    ### Import modules
    import numpy as np
    import calc_Utilities as UT
    import read_CTLNQ as CONT
    import read_ExpMonthly as NUDG
    import read_ShortCoupled as COUP
    import read_SIT as THICK
    import read_SIC as CONC
    import read_SIT_E3SM as E3SIT
    import read_SIC_E3SM as E3SIC
    import read_OldIceExperi as OLD
    import read_LongCoupled as LC
    
    if any([vari=='T700',vari=='T500']):
        varia = 'TEMP'
        level = 'profile'
    elif vari == 'U700':
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
    if period == 'NONE':
        latq = np.where((lat >= latpolar))[0]
        anomp = anom[:,:,latq,:]
        lat2p = lat2[latq,:]
        polarave = UT.calc_weightedAve(anomp,lat2p)
    else:
        latq = np.where((lat >= latpolar))[0]
        anomp = anom[:,latq,:]
        lat2p = lat2[latq,:]
        polarave = UT.calc_weightedAve(anomp,lat2p)
    
    print('\n========Calculated Polar Cap Average========\n')
    return polarave

### Test functions (do not use!)
#ave = PolarCap('SIC_Pd','SIC','surface',30,'NONE')

#import matplotlib.pyplot as plt
#import numpy as np
#plt.figure(figsize=(11,4))
#plt.title('Monthly SIC Anomalies')
#plt.plot(ave.ravel())
#plt.savefig('/home/zlabe/Desktop/' + 'monthly_SIC_anom.png',dpi=300)