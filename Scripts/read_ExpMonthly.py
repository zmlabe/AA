"""
Script reads in monthly data from AA/UTW experiments with 51 ensembles
 
Notes
-----
    Author : Zachary Labe
    Date   : 14 February 2020
    
Usage
-----
    [1] readExperi(varid,exp,timeperiod,level,cps)
"""

def readExperi(varid,exp,timeperiod,level,cps):
    """
    Function reads monthly data for AA/UTW experiments

    Parameters
    ----------
    varid : string
        variable name to read
    level : string
        Height of variable (surface or profile)
    exp : string
        Name of experiment to compare for nudging
    timeperiod : string
        Name of decadal nudging
    cps : string
        none or yes for stratosphere

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
    lat,lon,lev,var = readExperi(varid,exp,timperiod,level,cps)
    """
    print('\n>>>>>>>>>> Using readExperi function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Call files for directory (1-51 members)
    if cps == 'none':
        directorydata = '/seley/ypeings/simu/'
        experi = exp + '-' + timeperiod
    elif cps == 'yes':
        directorydata = '/seley/ypeings/simu/'
        experi = exp + '-' + timeperiod + '-cps'
    
    totaldirectory = directorydata + experi + '/monthly/'
    filename = totaldirectory + varid + '_1950-2000.nc'

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
        lev = data.variables['level'][:]
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'zonmean': # 3d variables (zonal mean!)
        varidz = varid + '_' + level
        filename = totaldirectory + varidz + '_1950-2000.nc'
        data = Dataset(filename,'r')
        lev = data.variables['level'][:]
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        varq = data.variables['%s' % varid][:].squeeze()
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
        
    print('Completed: Read members 1-51!')

    print('>>>>>>>>>> Completed: Finished readExperi function!')
    return lat,lon,lev,var

###############################################################################
###############################################################################
### Test functions - do not use!
#lat,lon,lev,var = readExperi('U10','UTW','2090','surface','yes')