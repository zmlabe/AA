"""
Script reads in monthly data from SC-WACCM4 short coupled experiments from PAMIP
 
Notes
-----
    Author : Zachary Labe
    Date   : 14 February 2020
    
Usage
-----
    [1] readCOUPs(varid,timeperiod,level)
"""

def readCOUPs(varid,timeperiod,level):
    """
    Function reads monthly data from short coupled experiments in PAMIP

    Parameters
    ----------
    varid : string
        variable name to read
    level : string
        Height of variable (surface or profile)
    timeperiod : string
        PAMIP-2.1 or PAMIP-2.2 or PAMIP-2.3

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
    lat,lon,lev,var = readCOUPs(varid,timperiod,level)
    """
    print('\n>>>>>>>>>> Using readCOUPs function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Call files for directory (1-100 members)
    if timeperiod == 'C_Fu':
        experi = 'PAMIP-2.3'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
        print('-----------USING SC-WACCM4 SHORT COUPLED EXPERIMENTS (Future)!-----------')
    elif timeperiod == 'C_Pd':
        experi = 'PAMIP-2.1'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
        print('-----------USING SC-WACCM4 SHORT COUPLED EXPERIMENTS (Present-Day)!-----------')
    elif timeperiod == 'C_Pi':
        experi = 'PAMIP-2.2'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
        print('-----------USING SC-WACCM4 SHORT COUPLED EXPERIMENTS (Pre-Industrial)!-----------')

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
        filename = totaldirectory + varidz + '_1900-2000.nc'
        data = Dataset(filename,'r')
        lev = data.variables['level'][:17]
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
        
    print('Completed: Read members 1-100!')

    print('>>>>>>>>>> Completed: Finished readCOUPs function!')
    return lat,lon,lev,var

###############################################################################
###############################################################################
### Test functions - do not use!
#import numpy as np
#import matplotlib.pyplot as plt
#vv = 'Z500'
#lat,lon,lev,var1 = readCOUPs(vv,'C_Fu','surface')
#lat,lon,lev,var = readCOUPs(vv,'C_Pd','surface')
#lat,lon,lev,var2 = readCOUPs(vv,'C_Pi','surface')