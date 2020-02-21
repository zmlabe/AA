"""
Script reads in monthly data from SC-WACCM4 SIC experiments from PAMIP
 
Notes
-----
    Author : Zachary Labe
    Date   : 18 February 2020
    
Usage
-----
    [1] readSIC(varid,timeperiod,level)
"""

def readSIC(varid,timeperiod,level):
    """
    Function reads monthly data from SIC experiments from PAMIP

    Parameters
    ----------
    varid : string
        variable name to read
    timeperiod : string
        Fu or Pd
    level : string
        Height of variable (surface or profile)

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
    lat,lon,lev,var = readSIC(varid,timperiod,level)
    """
    print('\n>>>>>>>>>> Using readSIC function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Call files for directory (1-301 members)
    if timeperiod == 'Fu':
        experi = 'PAMIP-1.6-QBO-300yr'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1700-2000.nc'
        print('-----------USING SC-WACCM4 CONCENTRATION EXPERIMENTS (Future)!-----------')
    elif timeperiod == 'Pd':
        experi = 'PAMIP-1.1-QBO-300yr'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1700-2000.nc'
        print('-----------USING SC-WACCM4 CONCENTRATION EXPERIMENTS (Present-Day)!-----------')
        
    ### Missing variables in other seley directory
    if any([varid=='THICK']):
        directorydata = '/seley/zlabe/simu/' # this is different
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1700-2000.nc'
        print('***ZLABE DIRECTORY***')

    ### Read in Data
    if level == 'surface': # 3d variables
        data = Dataset(filename,'r')
        lev = 'surface'
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'profile': # 4d variables
        data = Dataset(filename,'r')
        lev = data.variables['level'][:17] # goes to 10 hPa instead of 1 hPa
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:,:17,:,:] # goes to 10 hPa instead of 1 hPa
        data.close()
    elif level == 'zonmean': # 3d variables (zonal mean!)
        varidz = varid + '_' + level
        filename = totaldirectory + varidz + '_1700-2000.nc'
        data = Dataset(filename,'r')
        lev = data.variables['level'][:17] # goes to 10 hPa instead of 1 hPa
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        varq = data.variables['%s' % varid][:,:17,:,:].squeeze() # goes to 10 hPa instead of 1 hPa
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
        
    print('Completed: Read members 1-301!')

    print('>>>>>>>>>> Completed: Finished readSIC function!')
    return lat,lon,lev,var

###############################################################################
###############################################################################
### Test functions - do not use!
#lat,lon,lev,var1 = readSIC('Z500','Fu','surface')
#lat,lon,lev,var = readSIC('Z500','Pd','surface')