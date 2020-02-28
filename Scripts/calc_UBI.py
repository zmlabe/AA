def UBI(simu,period,varia):
    """
    Script calculates the Ural Blocking Index (SLP)
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
    
    ### Select variable
    level = 'surface'
    
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
    elif simu == 'coupled':
        lat,lon,lev,future = COUP.readCOUPs(varia,'C_Fu',level)
        lat,lon,lev,historical = COUP.readCOUPs(varia,'C_Pd',level)        
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
    
    ### Meshgrid for lat,lon
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate UBI
    lonq1 = np.where((lon >=0) & (lon <=90))[0]
    lonq2 = np.where((lon >= 330) & (lon <= 360))[0]
    lonq = np.append(lonq1,lonq2)
    latq = np.where((lat >=45) & (lat <=80))[0]
    anomlon = anom[:,:,lonq]
    anomu = anomlon[:,latq,:]
    lat2uq = lat2[latq,:]
    lat2u = lat2uq[:,lonq]
    ubi = UT.calc_weightedAve(anomu,lat2u)
    
    print('\n========Calculated Ural Blocking Index=======\n')
    return ubi

### Test functions (do not use!)
#ubi = UBI('AA-2090','DJF','Z500')