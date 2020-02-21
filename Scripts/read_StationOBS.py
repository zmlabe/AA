def readStationData(datasets,years):
    """ Read in station-based data sets for T2M
    """
    
    ### Import modules
    import numpy as np
    
    ### Define directories
    directorydata = '/home/zlabe/Documents/Projects/IceVarFigs/Data/'
    
    ### Set time
    yearall = np.arange(1900,2018+1,1)
    
    ### Read in data
    datat = np.empty((len(datasets),len(yearall)))
    for i in range(len(datasets)):
        datat[i] = np.genfromtxt(directorydata + '%s_Arctic_%s.txt' % (datasets[i],
                                 2018),delimiter=',',skip_header=1,
                                 unpack=True,usecols=[1])
    
    ### Mask correct time series
    yearq = np.where((yearall >= years.min()) & (yearall <= years.max()))[0]
    temp = datat[:,yearq]

    return temp

def calcStationEpoch(data,epoch):
    """
    Calculate epoch differences
    """
    
    ### Import modules
    import numpy as np
    
    diff = np.empty((len(data)))
    for i in range(len(data)):
        old = np.nanmean(data[i,:epoch])
        new = np.nanmean(data[i,-epoch:])
        diff[i] = new - old
        
    return diff

### Call functions (do not use!)
#import numpy as np
#import matplotlib.pyplot as plt
#datasets = ['GISTEMP','Berkeley']
#data = readStationData(datasets,np.arange(1979,2017+1,1))
#diff = calcStationEpoch(data,10)