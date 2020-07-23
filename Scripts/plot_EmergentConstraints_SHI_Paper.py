"""
Script plots relationship between vertical warming and the Siberian High for 
all model simulations

Notes
-----
    Author : Zachary Labe
    Date   : 25 March 2020
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
import read_AMIPAA as AMAA
import read_AMIP6 as AM

### Define directories
directoryfigure = '/home/zlabe/Desktop/AA/Emergent/Paper_v1/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Scatter of Warming-High (ALL)- %s----' % titletime)

### Add parameters
datareader = False
latpolar = 65.
variable = 'THICK'
period = 'DJF' 
level = 'surface'
runnames = [r'$\Delta$AA-2030',r'$\Delta$AA-2060',r'$\Delta$AA-2090',
            r'$\Delta$S-Coupled-Pd',r'$\Delta$S-Coupled-Pi',r'$\Delta$WACCM-SIT-Pd',
            r'$\Delta$WACCM-SIC-Pd',r'$\Delta$WACCM-SIC-Pi',r'$\Delta$NET']
runnamesdata = ['AA-2030','AA-2060','AA-2090',
                'coupled_Pd','coupled_Pi','SIT',
                'SIC_Pd','SIC_Pi','OLD']

runnames_E3SM = [r'$\Delta$E3SM-Pi',r'$\Delta$E3SM-Pd']
runnamesdata_E3SM = ['E3SIC_Pi','E3SIC_Pd']

runnames_AMIP = [r'AMIP-HL',r'AMIP']
runnamesdata_AMIP = ['AMIP-HL','AMIP']

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
    years = np.arange(1979,2016+1,1)
    late,lone,leve,thicke = REAN.readOBS('ERA5',variable,level,period)
    late,lone,leve,slpe = REAN.readOBS('ERA5','SLP',level,period)
    she = REAN.calcOBS_SHI(slpe,late,lone)
    vare= REAN.calcOBS_PolarCap(thicke,late,lone,latpolar)
    
    latrr,lonrr,levr,thickr = REAN.readOBS('NCEP1',variable,level,period)
    latr,lonr,levr,slpr = REAN.readOBS('NCEP1','SLP',level,period)
    shr = REAN.calcOBS_SHI(slpr,latr,lonr)
    varr= REAN.calcOBS_PolarCap(thickr,latrr,lonrr,latpolar)
    
    ###########################################################################
    ### Read in AMIP data    
    yearsAA = np.arange(1979,2016+1,1)
    lataa,lonaa,timeaa,levaa,t_AA = AMAA.readAMIPAA(variable,'AMIP-AA',level,False,True,period)
    lataa,lonaa,timeaa,levaa,slp_AA = AMAA.readAMIPAA('SLP','AMIP-AA',level,False,True,period)
    t_AAe = np.nanmean(t_AA,axis=0)
    slp_AAe = np.nanmean(slp_AA,axis=0)
    sh_AA = REAN.calcOBS_SHI(slp_AAe,lataa,lonaa)
    varr_AA = REAN.calcOBS_PolarCap(t_AAe,lataa,lonaa,latpolar)
    
    lathl,lonhl,timehl,levhl,t_HL = AMAA.readAMIPAA(variable,'AMIP-HL',level,False,True,period)
    lathl,lonhl,timehl,levhl,slp_HL = AMAA.readAMIPAA('SLP','AMIP-HL',level,False,True,period)
    t_HLe = np.nanmean(t_HL,axis=0)
    slp_HLe = np.nanmean(slp_HL,axis=0)
    sh_HL = REAN.calcOBS_SHI(slp_HLe,lathl,lonhl)
    varr_HL = REAN.calcOBS_PolarCap(t_HLe,lathl,lonhl,latpolar)
    
    latreg,lonreg,timereg,levreg,t_reg = AM.readAMIP6(variable,'AMIP',level,False,True,period)
    latreg,lonreg,timereg,levreg,slp_reg = AM.readAMIP6('SLP','AMIP',level,False,True,period)
    t_rege = np.nanmean(t_reg,axis=0)
    slp_rege = np.nanmean(slp_reg,axis=0)
    sh_reg = REAN.calcOBS_SHI(slp_rege,lathl,lonhl)
    varr_reg = REAN.calcOBS_PolarCap(t_rege,lathl,lonhl,latpolar)
 
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

oldthickaa = np.nanmean(varr_AA[:epochq])
newthickaa = np.nanmean(varr_AA[-epochq:])
diffaa = newthickaa - oldthickaa

oldthickhl = np.nanmean(varr_HL[:epochq])
newthickhl = np.nanmean(varr_HL[-epochq:])
diffhl = newthickhl - oldthickhl

oldthickreg = np.nanmean(varr_reg[:epochq])
newthickreg = np.nanmean(varr_reg[-epochq:])
diffreg = newthickreg - oldthickreg

### Calculate reanalysis epochs for Siberian High
oldshe = np.nanmean(she[:epochq])  # 1979-1988
newshe = np.nanmean(she[-epochq:]) # 2008-2017
diffshe = newshe - oldshe

oldshr = np.nanmean(shr[:epochq])  # 1979-1988
newshr = np.nanmean(shr[-epochq:]) # 2008-2017
diffshr = newshr - oldshr

oldshaa = np.nanmean(sh_AA[:epochq])
newshaa = np.nanmean(sh_AA[-epochq:])
diffshaa = newshaa - oldshaa

oldshhl = np.nanmean(sh_HL[:epochq])
newshhl = np.nanmean(sh_HL[-epochq:])
diffshhl = newshhl - oldshhl

oldshreg = np.nanmean(sh_reg[:epochq])
newshreg = np.nanmean(sh_reg[-epochq:])
diffshreg = newshreg - oldshreg

### Combine AMIP runs
POLAMIP_ALL = [diffhl,diffreg]
SHIAMIP_ALL = [diffshhl,diffshreg]

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
meanPOL_ALL1 = np.append(meanPOL,meanPOL_E3SM)
meanPOL_ALL = np.append(meanPOL_ALL1,POLAMIP_ALL)
meanSHI_ALL1 = np.append(meanSHI,meanSHI_E3SM)
meanSHI_ALL = np.append(meanSHI_ALL1,SHIAMIP_ALL)
    
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
    plt.axhspan(diffshe,diffshr,alpha=1,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k')
    
    color = cmocean.cm.phase(np.linspace(0.0,0.93,len(runnames)-1))
#    color = plt.cm.cubehelix(np.linspace(0.0,1,len(runnames)))
    for i,c in zip(range(len(runnames)-1),color):
        if i < 3:
            plt.scatter(meanPOL[i],meanSHI[i],color=c,s=52,marker='s',
                        label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                        edgecolor='k',linewidth=0.5)
        else:
            plt.scatter(meanPOL[i],meanSHI[i],color=c,s=52,
                        label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                        edgecolor='k',linewidth=0.5)

    plt.scatter(meanPOL[8],meanSHI[8],color='gold',s=52,
                label=r'\textbf{%s}' % runnames[8],zorder=11,clip_on=False,
                edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp_r(np.linspace(0.3,0.7,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=52,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    color = cmocean.cm.dense_r(np.linspace(0.2,0.9,len(runnames_AMIP)))
    for i,c in zip(range(len(runnames_AMIP)),color):
        plt.scatter(POLAMIP_ALL[i],SHIAMIP_ALL[i],color=c,s=52,
                    label=r'\textbf{%s}' % runnames_AMIP[i],zorder=11,
                    clip_on=False,marker='v',edgecolor='k',linewidth=0.5)   
    leg = plt.legend(shadow=False,fontsize=7,loc='upper center',
                     bbox_to_anchor=(0.5,1.16),fancybox=True,ncol=5,frameon=False,
                     handlelength=0,handletextpad=1)
    
    plt.xticks(np.arange(0,100,10),map(str,np.arange(0,100,10)),size=8)
    plt.yticks(np.arange(-5,10,0.5),map(str,np.arange(-5,10,0.5)),size=8)
    plt.xlim([0,90])
    plt.ylim([-0.5,4])
    
    plt.xlabel(r'\textbf{$\bf{\Delta}$1000-500 Thickness [m]}',
                         color='k',size=11,labelpad=5)
    plt.ylabel(r'\textbf{$\bf{\Delta}$Siberian High Index [hPa]}',
                         color='k',size=11,labelpad=5)
    plt.text(91,-0.32,r'\textbf{R$\bf{^{2}}$ = %s' % np.round(r_value**2,2),
            color='k',ha='right')
    plt.text(91,-0.5,r'\textbf{\textit{P}$\bf{<}$0.001}',
            color='k',ha='right')
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_v1_%s.png' % variable,
                dpi=300)
    
elif variable == 'T700':
    plt.axvspan(diffe,diffr,alpha=1,color='dimgrey',clip_on=False)
    plt.axhspan(diffshe,diffshr,alpha=1,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k',clip_on=False)
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.3,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    color = cmocean.cm.rain(np.linspace(0.25,0.8,len(runnames_AMIP)))
    for i,c in zip(range(len(runnames_AMIP)),color):
        plt.scatter(POLAMIP_ALL[i],SHIAMIP_ALL[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_AMIP[i],zorder=11,
                    clip_on=False,marker='v',edgecolor='k',linewidth=0.5)   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.7),fancybox=True,ncol=1,frameon=False,
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
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_v1_%s.png' % variable,
                dpi=300)
 
elif variable == 'T500':
    plt.axvspan(diffe,diffr,alpha=1,color='dimgrey',clip_on=False)
    plt.axhspan(diffshe,diffshr,alpha=1,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k',clip_on=False)
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.3,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    color = cmocean.cm.rain(np.linspace(0.25,0.8,len(runnames_AMIP)))
    for i,c in zip(range(len(runnames_AMIP)),color):
        plt.scatter(POLAMIP_ALL[i],SHIAMIP_ALL[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_AMIP[i],zorder=11,
                    clip_on=False,marker='v',edgecolor='k',linewidth=0.5)   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.7),fancybox=True,ncol=1,frameon=False,
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
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_v1_%s.png' % variable,
                dpi=300)
    
elif variable == 'T2M':
    plt.axvspan(diffobs.min(),diffobs.max(),alpha=1,color='k',
                clip_on=False,linewidth=0)
    plt.axvspan(diffe,diffr,alpha=1,color='dimgrey',clip_on=False,linewidth=0)
    plt.axhspan(diffshe,diffshr,alpha=1,color='dimgrey',clip_on=False,linewidth=0)
    plt.plot(xaxis,linetrend,linewidth=2,color='k')
    
    color = cmocean.cm.thermal(np.linspace(0.01,1,len(runnames)))
    for i,c in zip(range(len(runnames)),color):
        plt.scatter(meanPOL[i],meanSHI[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames[i],zorder=11,clip_on=False,
                    edgecolor='k',linewidth=0.5)
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(0.935,0.3),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
        
    color = cmocean.cm.amp(np.linspace(0.3,1.1,len(runnames_E3SM)))
    for i,c in zip(range(len(runnames_E3SM)),color):
        plt.scatter(meanPOL_E3SM[i],meanSHI_E3SM[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_E3SM[i],zorder=11,clip_on=False,marker='x')   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.8),fancybox=True,ncol=1,frameon=False,
                     handlelength=0,handletextpad=1)
    
    color = cmocean.cm.rain(np.linspace(0.25,0.8,len(runnames_AMIP)))
    for i,c in zip(range(len(runnames_AMIP)),color):
        plt.scatter(POLAMIP_ALL[i],SHIAMIP_ALL[i],color=c,s=42,
                    label=r'\textbf{%s}' % runnames_AMIP[i],zorder=11,
                    clip_on=False,marker='v',edgecolor='k',linewidth=0.5)   
    leg = plt.legend(shadow=False,fontsize=8,loc='upper center',
                     bbox_to_anchor=(1,0.7),fancybox=True,ncol=1,frameon=False,
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
    
    plt.savefig(directoryfigure + 'SHI_EmergentConstraints_v1_%s.png' % variable,
                dpi=300)

