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
import scipy.stats as sts
import calc_PolarCap as CAP

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/Seasons/Coupled/Forcings/'

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
runnames = [r'$\Delta$S-Coupled-Pd',r'$\Delta$L-Coupled-Pd',r'$\Delta$WACCM-SIC-Pd',
            r'$\Delta$E3SM-SIC-Pd']
runnamesdata = ['coupled_Pd','LONG','SIC_Pd','E3SIC_Pd']

###############################################################################
###############################################################################
###############################################################################
### Read in data
if datareader == True:
    ###########################################################################
    ### Read in model data for SC-WACCM4
    meanpol = np.empty((len(runnamesdata),12))
    for i in range(len(runnames)):
        polaraveq = CAP.PolarCap(runnamesdata[i],variable,level,latpolar,period)
        meanpol[i,:] = np.nanmean(polaraveq,axis=0)

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

color = cmocean.cm.thermal(np.linspace(0.3,1,len(runnames)))
for i,c in zip(range(len(runnames)),color):
    plt.plot(meanpol[i],color=c,marker='o',
                label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                linewidth=2)
leg = plt.legend(shadow=False,fontsize=8,loc='upper left',
                 bbox_to_anchor=(0,0.3),fancybox=True,ncol=1,frameon=False,
                 handlelength=1,handletextpad=1)

plt.xticks(np.arange(0,12,1),map(str,np.arange(0,12,1)),size=8)
plt.yticks(np.arange(-30,1,1),map(str,np.arange(-30,1,1)),size=8)
plt.xlim([0,11])
plt.ylim([-5,0])

plt.xlabel(r'\textbf{Months}',
                     color='k',size=11,labelpad=5)
plt.ylabel(r'\textbf{$\bf{\Delta}$SIC [\%]}',
                     color='k',size=11,labelpad=5)

plt.savefig(directoryfigure + 'AnomaliesSeasonalCycle_%s.png' % variable,
            dpi=300)