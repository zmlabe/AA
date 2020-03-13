"""
Functions read in monthly data from the 1 AMIP experiments (AMIP-AA). Data is 
available over the 1979-2016 period (38 years) and sorted by month (12). The 
AMIP simulations use SC-WACCM4 with historical forcings and RCP 4.5. 
The experiment with all forcings is called AMQS. Note that the first year 
(1978) is removed due to model spin-up.
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 March 2020
    
Usage
-----
    [1] readAMIPAA(variable,experiment,level,detrend,sliceeq)
"""

def readAMIPAA(variable,experiment,level,detrend,sliceeq):
    """
    Function reads monthly data from all 1 AMIP experiment

    Parameters
    ----------
    variable : string
        variable name to read
    experiment : string
        experiment name (AMIP-AA)
    level : string
        Height of variable (surface or profile)
    detrend : binary
        True/False whether to remove a linear trend at all grid points
    sliceeq : binary
        True/False whether to slice at the equator for only northern hemisphere
        
        
    Returns
    -------
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    time : 1d numpy array
        standard time (months since 1978-1-15, 00:00:00)
    lev : 1d numpy array
        levels (17)
    var : 5d numpy array or 6d numpy array 
        [ensemble,year,month,lat,lon] or [ensemble,year,month,level,lat,lon]

    Usage
    -----
    lat,lon,time,lev,var = readDataM(variable,experiment,level,detrend)
    """
    print('\n>>> Using readDataM function! \n')
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    import calc_Detrend as DT
    
    ### Declare knowns
    ensembles = 10
    months = 12
    years = np.arange(1979,2016+1,1)
    
    ### Directory for experiments (remote server - Seley)
    directorydata = '/seley/ypeings/simu/AMIP-AA/'
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Read in lat,lon,time from known file 
    if level == 'surface': # 3d variables
        dataq = Dataset(directorydata + 'AMIP-AA1/monthly/T2M_1978-2016.nc')
        time = dataq.variables['time'][12:]
        lev = 'surface'
        lat = dataq.variables['latitude'][:]
        lon = dataq.variables['longitude'][:]
        dataq.close()
        
    ###########################################################################
    ###########################################################################                 
        if sliceeq == False:
            ### Create empty variable
            varq = np.empty((ensembles,time.shape[0],
                             lat.shape[0],lon.shape[0]))
            varq[:,:,:,:] = np.nan ### fill with nans
    
        elif sliceeq == True:
            ### Slice for Northern Hemisphere
            latq = np.where(lat >= 0)[0]
            lat = lat[latq]
            ### Create empty variable
            varq = np.empty((ensembles,time.shape[0],
                     lat.shape[0],lon.shape[0]))
            varq[:,:,:,:] = np.nan ### fill with nans
            print('SLICE for Northern Hemisphere!')
        else:
            print(ValueError('Selected wrong slicing!'))
    
    ###########################################################################
    ###########################################################################
    elif level == 'profile': # 4d variables
        dataq = Dataset(directorydata + 'AMIP-AA1/monthly/TEMP_1978-2016.nc')
        time = dataq.variables['time'][12:]
        lev = dataq.variables['level'][:]
        lat = dataq.variables['latitude'][:]
        lon = dataq.variables['longitude'][:]
        dataq.close()
        
    ###########################################################################
    ###########################################################################
        if sliceeq == False:
            ### Create empty variable
            varq = np.empty((ensembles,time.shape[0],lev.shape[0],
                             lat.shape[0],lon.shape[0]))
            varq[:,:,:,:,:] = np.nan ### fill with nans
        elif sliceeq == True:
            ### Slice for Northern Hemisphere
            latq = np.where(lat >= 0)[0]
            lat = lat[latq]
            ### Create empty variable
            varq = np.empty((ensembles,time.shape[0],lev.shape[0],
                             lat.shape[0],lon.shape[0]))
            varq[:,:,:,:,:] = np.nan ### fill with nans
            print('SLICE for Northern Hemisphere!')
        else:
            print(ValueError('Selected wrong slicing!'))
    
    ###########################################################################
    ###########################################################################
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))
    
    ###########################################################################
    ###########################################################################
    ### Path name for file for each ensemble member
    for i in range(ensembles):
        filename = directorydata + '%s%s/' % (experiment,i+1) + \
                    'monthly/' + variable + '_1978-2016.nc'
                    
    ###########################################################################
    ###########################################################################
        ### Read in Data
        if sliceeq == False:
            if level == 'surface': # 3d variables
                data = Dataset(filename,'r')
                varq[i,:,:,:] = data.variables[variable][12:468,:,:]
                data.close()
                print('Completed: Read data %s%s- %s!' % (experiment,
                                                          i+1,variable))   
            elif level == 'profile': # 4d variables
                data = Dataset(filename,'r')
                varq[i,:,:,:,:] = data.variables[variable][12:468,:,:,:]
                data.close()
                print('Completed: Read data %s%s- %s!' % (experiment,
                                                          i+1,variable))
            else:
                print(ValueError('Selected wrong height - (surface or profile!)!'))
                
    ###########################################################################
    ###########################################################################
        elif sliceeq == True:
            if level == 'surface': # 3d variables
                data = Dataset(filename,'r')
                varq[i,:,:,:] = data.variables[variable][12:468,latq,:]
                data.close()
                print('Completed: Read data %s%s- %s!' % (experiment,
                                                          i+1,variable))      
            elif level == 'profile': # 4d variables
                data = Dataset(filename,'r')
                varq[i,:,:,:,:] = data.variables[variable][12:468,:,latq,:]
                data.close()
                print('Completed: Read data %s%s- %s!' % (experiment,
                                                          i+1,variable))
                
            else:
                print(ValueError('Selected wrong height - (surface or profile!)!'))
        
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Reshape to split years and months
    if level == 'surface': # 3d variables
        var = np.reshape(varq,(ensembles,varq.shape[1]//12,months,
                               lat.shape[0],lon.shape[0]))
    elif level == 'profile': # 4d variables
        var = np.reshape(varq,(ensembles,varq.shape[1]//12,months,lev.shape[0],
                      lat.shape[0],lon.shape[0]))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('\nCompleted: Reshaped %s array!' % (variable))
    
    ### Save computer memory
    del varq
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Convert units
    if variable in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')
    elif variable == 'SWE':
        var = var*1000. # Meters to Millimeters 
        print('Completed: Changed units (m to mm)!')

    ###########################################################################
    ###########################################################################
    ###########################################################################    
    ### Missing data (fill value to nans)
    var[np.where(var <= -8.99999987e+33)] = np.nan
    print('Completed: Filled missing data to nan!')
    
    ### Detrend data if turned on
    if detrend == True:
        var = DT.detrendData(var,level,'monthly')
    
    print('\n>>> Completed: Finished readAMIPAA function!')
    return lat,lon,time,lev,var
        
### Test functions -- no need to use    
#import numpy as np
#import matplotlib.pyplot as plt
#import scipy.stats as sts
#variable = 'SLP'
#experiment = 'AMIP-AA'
#level = 'surface'
#detrend = False
#sliceeq = False
#    
#lat,lon,time,lev,var = readAMIPAA(variable,experiment,level,detrend,sliceeq)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        