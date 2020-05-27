"""
Script reads in monthly data from ERA-Interim and performs various calculations
 
Notes
-----
    Author : Zachary Labe
    Date   :20 February 2020
    
Usage
-----
    [1] readOBS(experi,varid,level,period)
    [2] calcOBS_SHI(var,lat,lon)
    [3] calcOBS_UBI(var,lat,lon)
    [4] calcOBS_PolarCap(var,lat,lon,latpolar)
"""

def readOBS(experi,varid,level,period):
    """
    Function reads monthly data from ERA-Interim reanalysis

    Parameters
    ----------
    experi : string
        type of reanalysis 
    varid : string
        variable name to read
    level : string
        Height of variable (surface or profile)
    period : string
        Months to calculate over

    Returns
    -------
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    var : 4d numpy array or 5d numpy array 
        [year,month,lat,lon] or [year,month,level,lat,lon]

    Usage
    -----
    lat,lon,lev,var = readOBS(experi,varid,level,period)
    """
    print('\n>>>>>>>>>> Using readOBS function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    import calc_Utilities as UT
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Call files for directory 
    if experi == 'ERAI_Present': # (1979-2018)
        directorydata = '/seley/zlabe/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1979-2018.nc'
    elif experi == 'NCEP1': # (1948-2019)
        directorydata = '/seley/zlabe/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1948-2019.nc'       
    elif experi == 'ERA5': # (1979-2019)
        directorydata = '/seley/zlabe/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1979-2019.nc'      

    ### Read in Data
    if level == 'surface': # 3d variables
        data = Dataset(filename,'r')
        lev = 'surface'
        levq = 1
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'profile': # 4d variables
        data = Dataset(filename,'r')
        lev = data.variables['level'][:]
        levq = lev.shape[0]
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:],varid))

    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 3d variables
        var = np.reshape(varq,(varq.shape[0]//months,months,
                              int(lat.shape[0]),int(lon.shape[0])))
    elif level == 'profile': # 4d variables
        var = np.reshape(varq,(varq.shape[0]//months,months,int(lev.shape[0]),
                      int(lat.shape[0]),int(lon.shape[0])))
    elif level == 'zonmean': # 3d variables (zonal mean!)
        var = np.reshape(varq,(varq.shape[0]//months,months,int(lev.shape[0]),
                      int(lat.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped %s array!' % (varid))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')
    elif varid == 'SWE':
        var = var*1000. # Meters to Millimeters 
        print('Completed: Changed units (m to mm)!')
    if experi in ('ERAI_Present'):
        if varid == 'SLP':
            var = var/100 # Pa to hPa
            print('Completed: Changed units (Pa to hPa)!')
        elif varid == 'THICK':
            var = var/9.81 # m^2 s^-2 to m
            print('Completed: Changed units for geopotential!')
        elif any([varid=='Z1000',varid=='Z925',varid=='Z850',varid=='Z700',
                  varid=='Z500',varid=='Z300',varid=='Z250',varid=='Z100',
                  varid=='Z50',varid=='Z30',varid=='Z10',varid=='Z1',
                  varid=='GEOP']):
            var = var/9.81 # m^2 s^-2 to m
            print('Completed: Changed units for geopotential!')
    elif experi in ('ERA5'):
        if varid == 'SLP':
            var = var/100 # Pa to hPa
            print('Completed: Changed units (Pa to hPa)!')
        
    print('Completed: Read years 1979-2017!')
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Select same time period for reanalysis data (1979-2017)
    if experi == 'NCEP1':
        years = np.arange(1948,2019+1,1)
        yearq = np.where((years >= 1979) & (years <= 2019))[0]
        var = var[yearq]
    elif experi == 'ERAI_Present':
        years = np.arange(1979,2018+1,1)
        yearq = np.where((years >= 1979) & (years <= 2018))[0]
        var = var[yearq]
    elif experi == 'ERA5':
        years = np.arange(1979,2019+1,1)
        yearq = np.where((years >= 1979) & (years <= 2019))[0]
        var = var[yearq]
    
    ###########################################################################
    ###########################################################################
    ########################################################################### 
    ### Calculate over period
    if period == 'OND':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,-3:],axis=1)
    elif period == 'D':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,-1:],axis=1)
    elif period == 'DJF':
        print('Calculating over %s months!' % period)
        varm = UT.calcDecJanFeb(var,var,lat,lon,level,levq) 
    elif period == 'JFM':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,0:3],axis=1)
    elif period == 'JF':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,0:2],axis=1)
    elif period == 'FMA':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,1:4],axis=1)
    elif period == 'FM':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,1:3],axis=1)
    elif period == 'J':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,0:1],axis=1)
    elif period == 'F':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,1:2],axis=1)
    elif period == 'M':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,2:3],axis=1)
    elif period == 'MA':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,2:4],axis=1)
    elif period == 'JJA':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var[:,5:8],axis=1)
    elif period == 'annual':
        print('Calculating over %s months!' % period)
        varm = np.nanmean(var,axis=1)
    elif period == 'none':
        print('Calculating over %s months!' % period)
        varm = var
    else:
        print(ValueError('Selected wrong month period!'))

    print('>>>>>>>>>> Completed: Finished readOBS function!')
    return lat,lon,lev,varm

def calcOBS_SHI(var,lat,lon):
    """
    Script calculates the Siberian High Index for reanalysis
    """
    
    ### Import modules
    import numpy as np
    import calc_Utilities as UT
    
    ### Meshgrid for lat,lon
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate SHI
    lonq = np.where((lon >=80) & (lon <=120))[0]
    latq = np.where((lat >=40) & (lat <=65))[0]
    varlon = var[:,:,lonq]
    anoms = varlon[:,latq,:]
    lat2sq = lat2[latq,:]
    lat2s = lat2sq[:,lonq]
    varshi = UT.calc_weightedAve(anoms,lat2s)
    
    print('\n========Calculated *OBS* Siberian High Index========\n')
    return varshi

def calcOBS_UBI(var,lat,lon):
    """
    Script calculates the Ural Blocking Index for reanalysis
    """
    
    ### Import modules
    import numpy as np
    import calc_Utilities as UT
    
    ### Meshgrid for lat,lon
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate UBI
    lonq1 = np.where((lon >=0) & (lon <=90))[0]
    lonq2 = np.where((lon >= 330) & (lon <= 360))[0]
    lonq = np.append(lonq1,lonq2)
    latq = np.where((lat >=45) & (lat <=80))[0]
    varlon = var[:,:,lonq]
    varu = varlon[:,latq,:]
    lat2uq = lat2[latq,:]
    lat2u = lat2uq[:,lonq]
    varubi = UT.calc_weightedAve(varu,lat2u)
    
    print('\n========Calculated *OBS* Ural Blocking Index=======\n')
    return varubi

def calcOBS_PolarCap(var,lat,lon,latpolar):
    """
    Script calculates the polar cap average for reanalysis
    """
    ### Import modules
    import numpy as np
    import calc_Utilities as UT
    
    ### Meshgrid for lat,lon
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate SHI
    latq = np.where((lat >= latpolar))[0]
    if var.ndim == 3:
        varp = var[:,latq,:]
        lat2p = lat2[latq,:]
    elif var.ndim == 4:
        varp = var[:,:,latq,:]
        lat2p = lat2[latq,:]
    varave = UT.calc_weightedAve(varp,lat2p)
    
    print('\n========Calculated *OBS* Polar Cap Average========\n')
    return varave

### Test functions (do not use!)
import numpy as np
import matplotlib.pyplot as plt
varaa = 'THICK'
#lat,lon,lev,var = readOBS('ERA5',varaa,'surface','DJF')
#shi = calcOBS_SHI(var,lat,lon)
#ubi = calcOBS_UBI(var,lat,lon)
#varave = calcOBS_PolarCap(var,lat,lon,67.)

#lat,lon,lev,var1 = readOBS('NCEP1',varaa,'surface','DJF')
#shi1 = calcOBS_SHI(var1,lat,lon)
#ubi1 = calcOBS_UBI(var1,lat,lon)
#varave1 = calcOBS_PolarCap(var1,lat,lon,67.)

#lat,lon,lev,var2 = readOBS('ERAI_Present',varaa,'surface','DJF')
#shi2 = calcOBS_SHI(var2,lat,lon)
#ubi2 = calcOBS_UBI(var2,lat,lon)
#varave2 = calcOBS_PolarCap(var2,lat,lon,67.)