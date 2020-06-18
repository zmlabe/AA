#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:51:39 2020

@author: zlabe
"""

"""
Script plots sea ice concentration annual cycle for present-day PAMIP 
experiments

Notes
-----
    Author : Zachary Labe
    Date   : 7 April 2020
"""

### Import modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import calc_Utilities as UT
from netCDF4 import Dataset

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/Seasons/Coupled/Forcings/'
directorydata = '/home/zlabe/Documents/Research/AA/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting SIC Seasonal Cycle - %s----' % titletime)

### Add parameters
datareader = True
latpolar = 40.
variable = 'SIC'
level = 'surface'
period = 'NONE'
runnamesdata = ['AA-Control','AA-2030','AA-2060','AA-2090',
                'PAMIP-1.1','PAMIP-1.6']
months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV',
          'DEC']

###############################################################################
###############################################################################
###############################################################################
### Read in data
if datareader == True:
    ###########################################################################
    ### Read in data
    data = Dataset(directorydata + 'SIC_LENS_1979-2008.nc')
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    sicCONT = data.variables['SIC'][:]
    data.close()
    
    data = Dataset(directorydata + 'SIC_LENS_2030.nc')
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    sic30 = data.variables['SIC'][:]
    data.close()
    
    data = Dataset(directorydata + 'SIC_LENS_2060.nc')
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    sic60 = data.variables['SIC'][:]
    data.close()
    
    data = Dataset(directorydata + 'SIC_LENS_2090.nc')
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    sic90 = data.variables['SIC'][:]
    data.close()
    
    data = Dataset(directorydata + 'SIC_PAMIP-1.1.nc')
    latp = data.variables['lat'][:]
    lonp = data.variables['lon'][:]
    sic11 = data.variables['SIC'][:]
    data.close()
    
    data = Dataset(directorydata + 'SIC_PAMIP-1.6.nc')
    latp = data.variables['lat'][:]
    lonp = data.variables['lon'][:]
    sic16 = data.variables['SIC'][:]
    data.close()
    
    dataall = [sicCONT,sic30,sic60,sic90,sic11,sic16]
    lonall = [lon,lon,lon,lon,lonp,lonp]
    latall = [lat,lat,lat,lat,latp,latp]

### Slice at 30N
datan = []
latn = []
for i in range(len(dataall)):
    latq = np.where(latall[i] >= 30)[0]
    latqq = latall[i][latq]
    latn.append(latqq)
    
    dataallmask = dataall[i]
    dataallmask[np.where(dataallmask > 100)] = np.nan
    dataallmask[np.where(dataallmask < 1)] = np.nan
    newsic = dataallmask
    newsicn = newsic[:,latq,:]
    datan.append(newsicn)
    
### Average over Arctic
meann = []
lat2 = []
lon2 = []
for i in range(len(datan)):
    lon2n,lat2n = np.meshgrid(lonall[i],latn[i])
    meanq = UT.calc_weightedAve(datan[i],lat2n)
    meann.append(meanq)
    lat2.append(lat2n)
    lon2.append(lon2n)
    

################################################################################
################################################################################
################################################################################
### Create subplots of sea ice anomalies 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

color = cmocean.cm.thermal(np.linspace(0,1,len(runnamesdata)))
for i,c in zip(range(len(runnamesdata)),color):
    if i < 4:
        plt.plot(meann[i],color=c,marker='o',
                label=r'\textbf{%s}' % runnamesdata[i],zorder=11,clip_on=False,
                linewidth=2)
    else:
        plt.plot(meann[i],color=c,marker='s',
                label=r'\textbf{%s}' % runnamesdata[i],zorder=11,clip_on=False,
                linewidth=2)
leg = plt.legend(shadow=False,fontsize=8,loc='upper left',
                 bbox_to_anchor=(0,0.35),fancybox=True,ncol=1,frameon=False,
                 handlelength=1,handletextpad=1)

plt.xticks(np.arange(0,13,1),months,size=8)
plt.yticks(np.arange(0,101,10),map(str,np.arange(0,101,10)),size=8)
plt.xlim([0,11])
plt.ylim([0,90])

plt.xlabel(r'\textbf{Months}',
                     color='k',size=11,labelpad=5)
plt.ylabel(r'\textbf{SIC [\%]}',
                     color='k',size=11,labelpad=5)

plt.savefig(directoryfigure + 'SeasonalCycle_ALL.png',dpi=900)