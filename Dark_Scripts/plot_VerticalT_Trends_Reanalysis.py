#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 10:13:21 2020

@author: zlabe
"""

"""
Script plots vertical profile of temperature trends

Notes
-----
    Author : Zachary Labe
    Date   : 7 May 2020
"""

### Import modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import calc_Utilities as UT
import scipy.stats as sts
import read_OBS as REAN

### Define directories
directoryfigure = '/home/zlabe/Documents/Projects/AA/Dark_Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Reanalysis vertical T - %s----' % titletime)

### Add parameters
datareader = True
latpolar = 65.
epochq = 10
variable = 'TEMP'
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m"]
period = 'DJF' 
level = 'profile'
dataERA = 'ERA5'

###############################################################################
###############################################################################
###############################################################################
### Read in data
if datareader == True:
    ###########################################################################
    ### Read in reanalysis data
    years = np.arange(1979,2019+1,1)
    late,lone,leve,vare = REAN.readOBS(dataERA,variable,level,period) 
    
    varpole = np.nanmean(vare,axis=3)
    
###############################################################################
###############################################################################
###############################################################################
### Calculate trends
trend = np.empty((varpole.shape[1],varpole.shape[2]))
pruns = np.empty((varpole.shape[1],varpole.shape[2]))
x = np.arange(varpole.shape[0])
for i in range(varpole.shape[1]):
    for j in range(varpole.shape[2]):
         trend[i,j],intercept,r_value,pruns[i,j],std_err = \
         sts.linregress(x,varpole[:,i,j])
         
### Decadal trend
trendd = trend*10.
trendno = trendd.copy()

pruns[np.where(pruns < 0.05)] = 1.
pruns[np.where(pruns != 1.)] = 0
trendd = trendd * pruns
trenddd = trendd.copy()
trenddd[np.where(trenddd == 0.)] = np.nan

############################################################################
############################################################################
############################################################################
#### Plot temperature profile
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='darkgrey')
plt.rc('ytick',color='darkgrey')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='black')

fig = plt.figure()
ax1 = plt.subplot(111)

### Set limits for contours and colorbars
limit = np.arange(-1,1.01,0.025)
barlim = np.arange(-1,2,1)
zscale = np.array([1000,925,850,700,500,300,200])
timeqe,levqe = np.meshgrid(late,leve)

ax1.spines['top'].set_color('k')
ax1.spines['right'].set_color('k')
ax1.spines['bottom'].set_color('k')
ax1.spines['left'].set_color('k')
ax1.spines['left'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='k')
ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='k')    
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

cs = plt.contourf(timeqe,levqe,trenddd,limit,extend='both')
cs1 = plt.contourf(timeqe,levqe,trendno,limit,extend='both',
                  alpha=0.3,antialiased=True)
#cs1 = plt.contour(timeqe,levqe,trendd,np.arange(-5,5.5,0.2),
#                   colors='dimgrey')
#cs2 = plt.contourf(timeqe,levqe,prunse,colors='None',
#               hatches=['//////'],linewidths=0.4)

plt.yscale('log',nonposy='clip')
plt.ylim([1000,200])
plt.xticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=6)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
plt.xlim([45,90])
plt.minorticks_off()
plt.gca().invert_yaxis()

plt.axhline(leve[-1],linewidth=2,linestyle='--',dashes=(1,0.5),
            color='w')
plt.axhline(leve[-2],linewidth=2,linestyle='--',dashes=(1,0.5),
            color='k')
plt.axhline(leve[-3],linewidth=2,linestyle='--',dashes=(1,0.5),
            color='k')
plt.axhline(leve[-4],linewidth=2,linestyle='--',dashes=(1,0.5),
            color='k')
plt.axhline(leve[-6],linewidth=2,linestyle='--',dashes=(1,0.5),
            color='k')
plt.axhline(leve[-8],linewidth=2,linestyle='--',dashes=(1,0.5),
            color='k')
           
cs.set_cmap(cmocean.cm.balance)
cs1.set_cmap(cmocean.cm.balance)

plt.tight_layout()
cbar_ax = fig.add_axes([0.312,0.05,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001,labelsize=6)
cbar.outline.set_edgecolor('darkgrey')

plt.subplots_adjust(bottom=0.2)

plt.savefig(directoryfigure + 'VerticalT_Trends_Reanalysis_ERA5.png',dpi=600)
