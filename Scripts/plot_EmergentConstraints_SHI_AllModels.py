"""
Script plots relationship between vertical warming and the Siberian High for 
all model simulations

Notes
-----
    Author : Zachary Labe
    Date   : 27 February 2020
"""

### Import modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import calc_Utilities as UT
import scipy.stats as sts
import calc_SHI as SH
import calc_PolarCap as CAP
import read_OBS as REAN
import read_StationOBS as SOBS

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/Emergent/SHI/ALL/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Scatter of Warming-High (ALL)- %s----' % titletime)

### Add parameters
datareader = True
latpolar = 65.
variable = 'T500'
period = 'DJF' 
level = 'surface'
runnames = [r'AA-2030',r'AA-2060',r'AA-2090',
            r'2.3--2.1',r'$\Delta$SIT-Pd',r'$\Delta$SIC-Pi',r'$\Delta$SIC-Pd',r'$\Delta$NET']
runnamesdata = ['AA-2030','AA-2060','AA-2090','coupled','SIT','SIC_Pi','SIC_Pd','OLD']

runnames_E3SM = ['E3SM-SIT-Pd','E3SM-Pi','E3SM-Pd']
runnamesdata_E3SM = ['E3SIT','E3SIC_Pi','E3SIC_Pd']

###############################################################################
###############################################################################
###############################################################################
### Read in data
if datareader == True:
    ###########################################################################
    ### Read in model data for SC-WACCM4
    polarave = []
    high = []
    for i in range(len(runnames)):
        polaraveq = CAP.PolarCap(runnamesdata[i],variable,level,latpolar,period)
        highq = SH.SHI(runnamesdata[i],period)
        polarave.append(polaraveq)
        high.append(highq)
    ###########################################################################
    ### Read in model data for E3SM
    polarave_E3SM = []
    high_E3SM = []
    for i in range(len(runnames_E3SM)):
        polaraveq_E3SM = CAP.PolarCap(runnamesdata_E3SM[i],variable,level,latpolar,period)
        highq_E3SM = SH.SHI(runnamesdata_E3SM[i],period)
        polarave_E3SM.append(polaraveq_E3SM)
        high_E3SM.append(highq_E3SM)
    ###########################################################################
    ### Read in reanalysis data
    years = np.arange(1979,2017+1,1)
    late,lone,leve,thicke = REAN.readOBS('ERAI_Present',variable,level,period)
    late,lone,leve,slpe = REAN.readOBS('ERAI_Present','SLP',level,period)
    she = REAN.calcOBS_SHI(slpe,late,lone)
    vare= REAN.calcOBS_PolarCap(thicke,late,lone,latpolar)
    
    latr,lonr,levr,thickr = REAN.readOBS('NCEP1',variable,level,period)
    latr,lonr,levr,slpr = REAN.readOBS('NCEP1','SLP',level,period)
    shr = REAN.calcOBS_SHI(slpr,latr,lonr)
    varr= REAN.calcOBS_PolarCap(thickr,latr,lonr,latpolar)
 
###############################################################################
###############################################################################
###############################################################################
### Calculate ensemble means for WACCM4
meanPOL = np.empty((len(polarave)))
meanSHI = np.empty((len(high)))
for i in range(len(runnames)):
    meanPOL[i] = np.nanmean(polarave[i])
    meanSHI[i] = np.nanmean(high[i])
    
### Calculate ensemble means for E3SM
meanPOL_E3SM = np.empty((len(polarave_E3SM)))
meanSHI_E3SM = np.empty((len(high_E3SM)))
for i in range(len(runnames_E3SM)):
    meanPOL_E3SM[i] = np.nanmean(polarave_E3SM[i])
    meanSHI_E3SM[i] = np.nanmean(high_E3SM[i])

###############################################################################
###############################################################################
###############################################################################
### Calculate reanalysis epochs for temperature
epochq = 10
oldthicke = np.nanmean(vare[:epochq])  # 1979-1988
newthicke = np.nanmean(vare[-epochq:]) # 2008-2017
diffe = newthicke - oldthicke

oldthickr = np.nanmean(varr[:epochq])  # 1979-1988
newthickr = np.nanmean(varr[-epochq:]) # 2008-2017
diffr = newthickr - oldthickr

### Calculate reanalysis epochs for Siberian High
oldshe = np.nanmean(she[:epochq])  # 1979-1988
newshe = np.nanmean(she[-epochq:]) # 2008-2017
diffshe = newshe - oldshe

oldshr = np.nanmean(shr[:epochq])  # 1979-1988
newshr = np.nanmean(shr[-epochq:]) # 2008-2017
diffshr = newshr - oldshr

### Calculate T2M station-based data sets
if variable == 'T2M':
    dataobs = SOBS.readStationData(['GISTEMP','Berkeley'],years)
    diffobs = SOBS.calcStationEpoch(dataobs,epochq)

################################################################################
################################################################################
################################################################################
#### Calculate statistics
if variable == 'THICK':
    xaxis = np.arange(0,100,10)
elif variable == 'T500':
    xaxis = np.arange(0,2.1,0.1)
elif variable == 'T700':
    xaxis = np.arange(0,4.2,0.1)
elif variable == 'T2M':
    xaxis = np.arange(0,10,0.1)

### Combine all model data
meanPOL_ALL = np.append(meanPOL,meanPOL_E3SM)
meanSHI_ALL = np.append(meanSHI,meanSHI_E3SM)
    
### Calculate regression
slope,intercept,r_value,p_value,std_err = sts.linregress(meanPOL_ALL,meanSHI_ALL)
linetrend = slope*xaxis + intercept

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
        
fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

if variable == 'THICK':
    plt.axvspan(diffe,diffr,alpha=1,color='dimgrey',clip_on=False)
    plt.axhspan(diffshe,diffshr,alpha=0.6,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k')
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.15,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    plt.xticks(np.arange(0,100,10),map(str,np.arange(0,100,10)),size=8)
    plt.yticks(np.arange(0,10,0.5),map(str,np.arange(0,10,0.5)),size=8)
    plt.xlim([0,90])
    plt.ylim([0,4])
    
    plt.xlabel(r'\textbf{$\bf{\Delta}$1000-500 Thickness [m]}',
                         color='k',size=11,labelpad=5)
    plt.ylabel(r'\textbf{$\bf{\Delta}$Siberian High Index [hPa]}',
                         color='k',size=11,labelpad=5)
    plt.text(91,3.9,r'\textbf{R$\bf{^{2}}$=%s' % np.round(r_value**2,2),
            color='k',ha='right')
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_PAMIP-Nudge_%s.png' % variable,
                dpi=300)
    
elif variable == 'T700':
    plt.axvspan(diffe,diffr,alpha=1,color='dimgrey',clip_on=False)
    plt.axhspan(diffshe,diffshr,alpha=0.6,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k',clip_on=False)
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.15,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    plt.xticks(np.arange(0,7.01,0.5),map(str,np.arange(0,7.01,0.5)),size=8)
    plt.yticks(np.arange(0,10,0.5),map(str,np.arange(0,10,0.5)),size=8)
    plt.xlim([0,4])
    plt.ylim([0,4])
    
    plt.xlabel(r'\textbf{$\bf{\Delta}$T700 [$\bf{^{\circ}}$C]}',
                         color='k',size=11,labelpad=5)
    plt.ylabel(r'\textbf{$\bf{\Delta}$Siberian High Index [hPa]}',
                         color='k',size=11,labelpad=5)
    
    plt.text(0,3.9,r'\textbf{R$\bf{^{2}}$=%s' % np.round(r_value**2,2),
            color='k')
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_PAMIP-Nudge_%s.png' % variable,
                dpi=300)
 
elif variable == 'T500':
    plt.axvspan(diffe,diffr,alpha=1,color='dimgrey',clip_on=False)
    plt.axhspan(diffshe,diffshr,alpha=0.6,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k',clip_on=False)
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.15,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    plt.xticks(np.arange(0,2.01,0.25),map(str,np.arange(0,2.01,0.25)),size=8)
    plt.yticks(np.arange(0,10,0.5),map(str,np.arange(0,10,0.5)),size=8)
    plt.xlim([0,2])
    plt.ylim([0,4])
    
    plt.xlabel(r'\textbf{$\bf{\Delta}$T500 [$\bf{^{\circ}}$C]}',
                         color='k',size=11,labelpad=5)
    plt.ylabel(r'\textbf{$\bf{\Delta}$Siberian High Index [hPa]}',
                         color='k',size=11,labelpad=5)
    
    plt.text(0,3.9,r'\textbf{R$\bf{^{2}}$=%s' % np.round(r_value**2,2),
            color='k')
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_PAMIP-Nudge_%s.png' % variable,
                dpi=300)
    
elif variable == 'T2M':
    plt.axvspan(diffobs.min(),diffobs.max(),alpha=0.6,color='darkgrey',
                clip_on=False,linewidth=0)
    plt.axvspan(diffe,diffr,alpha=0.6,color='dimgrey',clip_on=False,linewidth=0)
    plt.axhspan(diffshe,diffshr,alpha=0.6,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k')
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.15,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    plt.xticks(np.arange(0,16,2),map(str,np.arange(0,16,2)),size=8)
    plt.yticks(np.arange(-10,10,0.5),map(str,np.arange(-10,10,0.5)),size=8)
    plt.xlim([0,10])
    plt.ylim([-1,4])
    
    plt.xlabel(r'\textbf{$\bf{\Delta}$T2M [$\bf{^{\circ}}$C]}',
                         color='k',size=11,labelpad=5)
    plt.ylabel(r'\textbf{$\bf{\Delta}$Siberian High Index [hPa]}',
                         color='k',size=11,labelpad=5)
    
    plt.text(0,3.9,r'\textbf{R$\bf{^{2}}$=%s' % np.round(r_value**2,2),
            color='k')
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_ALL_%s.png' % variable,
                dpi=300)

