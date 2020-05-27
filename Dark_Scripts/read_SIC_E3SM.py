"""
Script reads in monthly data from E3SM SIC experiments from PAMIP
 
Notes
-----
    Author : Zachary Labe
    Date   : 27 February 2020
    
Usage
-----
    [1] readE3SM_SIC(varid,timeperiod,level)
"""

def readE3SM_SIC(varid,timeperiod,level):
    """
    Function reads monthly data from SIC experiments from PAMIP (E3SM)

    Parameters
    ----------
    varid : string
        variable name to read
    level : string
        Height of variable (surface or profile)
    timeperiod : string
        ESIC_Fu or ESIC_Pd or ESIC_Pi

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
    lat,lon,lev,var = readE3SM_SIC(varid,timperiod,level)
    """
    print('\n>>>>>>>>>> Using readE3SM_SIC function!')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Call files for directory (1-101 members)
    if timeperiod == 'ESIC_Fu':
        experi = 'PAMIP-1.6-E3SM'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
        print('-----------USING DOE E3SM *CONC* EXPERIMENTS!-----------')
    elif timeperiod == 'ESIC_Pd':
        experi = 'PAMIP-1.1-E3SM'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
        print('-----------USING DOE E3SM *CONC* EXPERIMENTS!-----------')
    elif timeperiod == 'ESIC_Pi':
        experi = 'PAMIP-1.5-E3SM'
        directorydata = '/seley/ypeings/simu/'
        totaldirectory = directorydata + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
        print('-----------USING DOE E3SM *CONC* EXPERIMENTS!-----------')
    else:
        print(ValueError('Selected wrong time period (ESIC_Fu,ESIC_Pd,ESIC_Pi!')) 

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
        
    print('Completed: Read members 1-101!')

    print('>>>>>>>>>> Completed: Finished readE3SM_SIC function!')
    return lat,lon,lev,var

###############################################################################
###############################################################################
### Test functions - do not use!
#lat,lon,lev,var2 = readE3SM_SIC('Z500','ESIC_Pi','surface')
#lat,lon,lev,var1 = readE3SM_SIC('Z500','ESIC_Fu','surface')
#lat,lon,lev,var = readE3SM_SIC('THICK','ESIC_Pd','surface')