# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 15:11:28 2021
Plot the modeling result for fig 2. Need to run Canopy16017034 first to generate required files first
@author: Hongyu Meng
"""
import matplotlib.pyplot as plt

from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np




MyColors = np. array([
    [80 ,      255,      80] ,   #Green       Ir
    [255,      255,      80] ,   # Yellow      ramp
    [255,      80 ,      80] ,   # Red         ib
    [80 ,      80 ,      255],   # Blue        adap
    
    [167 ,     255,      80] ,   # grass-gree? ir/ramp
    [255,      200,      80 ],  # orange      ir/ib  
    [80 ,      255,      255],   # cyan        ir/adap  
    [255,      80 ,      255],   # purple      ib/adap
    [167,      167,      167],   # grey      other
    [255,      255 ,     255],   # write      all
    ]);  

MyColors=MyColors*180/255

MyColors=MyColors/255

ColorIndex=0






# Plot the single traces. 
checkregion=b2.array([-10, 100])
CheckInd=160  # 49 for 200 pA


# For better plotting, I need to plot the voltage to the -20 mV

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))



spike_train=spike_trains[CheckInd]
Ind_spike=(spike_train/b2.defaultclock.dt).astype(int)

vcopy=statemon.v[CheckInd]/b2.mV
vcopy[Ind_spike]=-20
 
#plt.title(r'phase diagram')

#ax.plot(statemon.t/b2.ms, statemon.v[CheckInd]/b2.mV, '-k',linewidth=1)
ax.plot(statemon.t/b2.ms-1000, vcopy, '-k',linewidth=1)   # The injection starts at 0
ax.set_xlim(checkregion)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([-70, -20])

plt.show()


    
fig.savefig(('Volt%2.0f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')



CheckInd=161 
Ind_start=1;
#Ind_end=np.int(simulation_time/(b2.defaultclock.dt));   # 10000 for 1s
Ind_end=np.int(2000*b2.ms/(b2.defaultclock.dt));  
volt_trace=statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV
Dif_trace= (statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV-statemon.v[CheckInd,Ind_start-1:Ind_end-1]/b2.mV)/(b2.defaultclock.dt/b2.ms)
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

ax.plot(volt_trace,Dif_trace, '-k',linewidth=1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

ax.set_xlim([-50,-30])
ax.set_ylim([-2, 20])
fig.savefig(('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')


########## Plot the fI curve
xlimrange=np.array([150, 300])

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
b2.plot(Inj_vec/b2.pA, spikemon.count/(simulation_time),color='black')
ax.set_xlim([150,300])
ax.set_ylim([0, 60] )
plt.xticks(np.arange(150, 300+1, 50.0))


#b2.xlabel('Injection (nA)')
#b2.ylabel('Firing rate (sp/s)')

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_ylim([0, 1] )

fig.savefig('Ficurve'+'.svg',bbox_inches='tight')

######  Plot different fI curves

'''
    spike_train_all=spike_trains[CheckInd]/b2.ms
    spike_train=spike_train_all[spike_train_all<500+1000]

    if  spike_train.size>=2:
        # Non-linear Fit
        ISI_vec=spike_train[1:]-spike_train[:-1]
        FI_vec=1000/(ISI_vec)
        ax.plot(spike_train[:-1]-500, FI_vec, '.',markersize=3,color=MyColors[ColorIndex,:] )   # The injection starts at 0

        
    if  spike_train.size>=6:
        A, K, C = fit_exp_nonlinear(spike_train[:-1]-500, 1000/(ISI_vec), [-FI_vec[-1]/6, -0.001, FI_vec[-1]]) 
        fit_y = model_func(spike_train[:-1]-500, A, K, C)
        ax.plot(spike_train[:-1]-500, fit_y, '-',linewidth=1,color=MyColors[ColorIndex,:]) 
    ColorIndex+=1    
ax.set_xlim([0, 1000])
#ax.set_xlim([0, 2000])

ax.set_ylim([0, 30])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(np.arange(0, 1000+1, 500.0))

fig.savefig('DifRateTrace'+'.svg',bbox_inches='tight')

'''


spike_trains = spikemon.spike_trains()

ColorIndex=0

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

for i_ind in range(5):
    CheckInd=i_ind*10+152
    
    spike_train=spike_trains[CheckInd]/b2.ms
    ISI_vec=spike_train[1:]-spike_train[:-1]

    ax.plot(spike_train[:-1]-1000, 1000/(ISI_vec), '.k',markersize=3,color=MyColors[ColorIndex,:] )   # The injection starts at 0
    ax.plot(spike_train[:-1]-1000, 1000/(ISI_vec), linewidth=1,color=MyColors[ColorIndex,:] )   # The injection starts at 0
    ColorIndex+=1
    
ax.set_xlim([0,100])
ax.set_ylim([0, 350])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(np.arange(0, 100+1, 20.0))

fig.savefig('DifRateTrace'+'.svg',bbox_inches='tight')





############# Plot std over injection current

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,Std_vec,color='black')
ax.set_xlim(xlimrange )
plt.xticks(np.arange(150, 300+1, 50.0))

#plt.xticks(np.arange(50, 250+1, 100.0))

plt.yticks(np.arange(0, 1.0+0.1, 0.5))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([0, 1] )

fig.savefig('VoltStd'+'.svg',bbox_inches='tight')
