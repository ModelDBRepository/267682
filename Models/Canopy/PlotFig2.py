# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 15:11:28 2021
Plot the modeling result for fig 2. Need to run Canopy16017034 first to generate required files first
@author: Hongyu Meng
"""
import matplotlib.pyplot as plt
import os
import brian2 as b2

from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

TempPath='Fig2Plots'
if os.path.isdir(TempPath)==False:
    os.makedirs(TempPath)
# Plot the single traces. 
    
    
checkregion=b2.array([-100, 2000])

Inj_vec= b2.linspace(1,Num_neu,Num_neu)/Num_neu *(Inj_max-Inj_min)+Inj_min  

# For better plotting, I need to plot the voltage to the -20 mV

CheckInd=79  # 30 

###########   Volt Trace
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   

xlimrange=np.array([300 ,400])

spike_train=spike_trains[CheckInd]
Ind_spike=(spike_train/b2.defaultclock.dt).astype(int)

vcopy=statemon.v[CheckInd]/b2.mV
vcopy[Ind_spike]=-20
 
#plt.title(r'phase diagram')

#ax.plot(statemon.t/b2.ms, statemon.v[CheckInd]/b2.mV, '-k',linewidth=1)
ax.plot(statemon.t/b2.ms-500, vcopy, '-k',linewidth=1)   # The injection starts at 0
ax.set_xlim(checkregion)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.show()

fig.savefig(TempPath+'/'+('Volt%.3f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')



###########   dVolt/dt Trace
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   
Ind_start=np.int((simulation_time-1.0*b2.second)/(b2.defaultclock.dt));
Ind_end=np.int(simulation_time/(b2.defaultclock.dt));   # 10000 for 1s

volt_trace=statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV
Dif_trace= (statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV-statemon.v[CheckInd,Ind_start-1:Ind_end-1]/b2.mV)/(b2.defaultclock.dt/b2.ms)

ax.plot(volt_trace,Dif_trace,'-k',linewidth=1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim([-55,-30])
ax.set_ylim([-2, 10])

ax.set_xlim([-60,-20])
ax.set_ylim([-300, 300])
b2.savefig(TempPath+'/'+('PhasePlanENlargedInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')

#b2.savefig(TempPath+'/'+('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg')

###########   a, b Trace
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   

ax.plot(statemon.t/b2.ms-500, statemon.a[CheckInd],'k',label='a',linewidth=1)
ax.plot(statemon.t/b2.ms-500, statemon.b[CheckInd],'r',label='b',linewidth=1)


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim(checkregion)
ax.set_ylim([0, 1])

b2.savefig(TempPath+'/'+('abInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')

#b2.savefig(TempPath+'/'+('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg')


#####################################
#### Plotting the max/min and CV_ISI in a twin axis.

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

ax.plot(Inj_vec/b2.pA,1-1/MaxMinratio_vec,color='black')
#ax.set(ylabel='CV(ISI)')
ax.set_ylim([0, 1] )

ax.spines["top"].set_visible(False)
ax.set_xlim(xlimrange )
plt.yticks(np.arange(0, 1.0+0.1, 0.5))

ax2=ax.twinx()
ax2.spines["top"].set_visible(False)

ax2.plot(Inj_vec/b2.pA,CVISI_vec,color='red')
#ax2.set(ylabel='Max/Min')
ax2.set_ylim([0, 1] )

ax2.yaxis.label.set_color('red')

ax2.tick_params(axis='y', colors='red')
ax2.spines['right'].set_color('red')

ax2.set_xlim(xlimrange )
plt.xticks(np.arange(300, 400+1, 50.0))
plt.yticks(np.arange(0, 1.0+0.1, 0.5))

fig.savefig(TempPath+'/'+'CVISIMaxmin'+'.svg',bbox_inches='tight')

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

ax.plot(Inj_vec/b2.pA,CVISI_vec,color='black')
#ax.set(ylabel='CV(ISI)')
ax.set_ylim([0, 1] )

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim(xlimrange )
plt.xticks(np.arange(300, 400+1, 50.0))
plt.yticks(np.arange(0, 1.0+0.1, 0.5))


fig.savefig(TempPath+'/'+'CVISI'+'.svg',bbox_inches='tight')

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

ax.plot(Inj_vec/b2.pA,1-1/MaxMinratio_vec,color='black')
#ax.set(ylabel='CV(ISI)')
ax.set_ylim([0, 1] )

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim(xlimrange )
plt.xticks(np.arange(300, 400+1, 50.0))
plt.yticks(np.arange(0, 1.0+0.1, 0.5))


fig.savefig(TempPath+'/'+'ISIratio'+'.svg',bbox_inches='tight')

# Plotting the a_inf; b_inf for the fig 2.
Num_v=100
v_max=-20*b2.mV; v_min=-60*b2.mV

v_vec=np.linspace(1,Num_v,Num_v)/Num_v *(v_max-v_min)+v_min
#mid_a=-37;sig_a=10
#mid_b=-52;sig_b=-6

ainf_vec=1/(1+np.exp(-(v_vec- mid_a) /sig_a      )   )
binf_vec=1/(1+np.exp(-(v_vec- mid_b) /sig_b      )   )


fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

ax.plot(v_vec*1000,ainf_vec,color='black')
ax.plot(v_vec*1000,binf_vec,color='red')
ax.set_ylim([0, 1] )
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.savefig(TempPath+'/'+'abinf'+'.svg',bbox_inches='tight')


fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,Std_vec,color='black')
ax.set_xlim(xlimrange )
plt.xticks(np.arange(300, 400+1, 50.0))
plt.yticks(np.arange(0, 1.0+0.1, 0.5))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([0, 1] )

fig.savefig(TempPath+'/'+'VoltStd'+'.svg',bbox_inches='tight')

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
b2.plot(Inj_vec/b2.pA, spikemon.count/(simulation_time),color='black')
ax.set_xlim(xlimrange )
ax.set_ylim([0, 60] )
plt.xticks(np.arange(300, 400+1, 50.0))


#b2.xlabel('Injection (nA)')
#b2.ylabel('Firing rate (sp/s)')

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_ylim([0, 1] )

fig.savefig(TempPath+'/'+'Ficurve'+'.svg',bbox_inches='tight')