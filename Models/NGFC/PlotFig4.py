# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 15:11:28 2021
Plot the modeling result for fig 2. Need to run Canopy16017034 first to generate required files first
@author: Hongyu Meng
"""
import matplotlib.pyplot as plt

from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import brian2 as b2
import scipy as sp
import scipy.optimize


def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def fit_exp_nonlinear(t, y,p0=[-5, 0.01, 25]):
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, t, y, p0, maxfev=10000)
    A, K, C = opt_parms
    return A, K, C


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
checkregion=b2.array([-100, 1000])
CheckInd=47  # 49 for 200 pA

# For better plotting, I need to plot the voltage to the -20 mV

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))



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


    
fig.savefig(('Volt%2.0f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')



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

#ax.set_xlim([-60,-20])
#ax.set_ylim([-300, 300])
#b2.savefig(('PhasePlanENlargedInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')
b2.savefig(('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')




############# Plot std over injection current
xlimrange=np.array([180, 210])

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,Std_vec,color='black')
ax.set_xlim([180,190] )
plt.xticks(np.arange(180, 190+0.1, 5.0))



plt.yticks(np.arange(0, 1.0+0.1, 0.5))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([0, 1] )

fig.savefig('VoltStd'+'.svg',bbox_inches='tight')


############# Plot 1st spike over injection current

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,First_spike,color='black')
ax.set_xlim(xlimrange )
ax.set_ylim([0, 800] )

plt.xticks(np.arange(180, 210+1, 10.0))
plt.yticks(np.arange(0, 800.0+0.1, 400.0))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_ylim([0, 1] )

fig.savefig('Time1stSpike'+'.svg',bbox_inches='tight')




######  Plot different fI curves

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ColorIndex=0
for i_ind in range(6):
    if i_ind==0:
        CheckInd=37
    else:
        CheckInd=35+i_ind*4
    
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

ax.set_ylim([0, 50])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(np.arange(0, 1000+1, 500.0))

fig.savefig('DifRateTrace'+'.svg',bbox_inches='tight')

FitInitRate=np.full(100,np.nan)     
FitLastRate=np.full(100,np.nan)     


 
for i_ind in range(100):
    
    spike_train_all=spike_trains[i_ind]/b2.ms
    spike_train=spike_train_all[spike_train_all<=500+1000]
#    ISI_vec=spike_train[1:]-spike_train[:-1]

    if  spike_train.size>=4:
        ISI_vec=spike_train[1:]-spike_train[:-1]
        FI_vec=1000/(ISI_vec)
#        ax.plot(spike_train[:-1]-500, 1000/(ISI_vec), '.k',markersize=1)   # The injection starts at 0
        # Non-linear Fit
        A, K, C = fit_exp_nonlinear(spike_train[:-1]-500, 1000/(ISI_vec), [-4, -0.00, FI_vec[-1]])   # last is p0. 
        fit_y = model_func(spike_train[:-1]-500, A, K, C)
#        ax.plot(spike_train[:-1]-500, fit_y, '-',linewidth=1) 
        FitInitRate[i_ind] =fit_y[1]
        FitLastRate[i_ind] =fit_y[-1]
#    ax.plot(spike_train[:-1]-500, 1000/(ISI_vec), linewidth=1)   # The injection starts at 0



############# Plot 1st IF and last IF over injection current

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,FitInitRate,color='black')
ax.plot(Inj_vec/b2.pA,FitLastRate,color='red')

ax.set_xlim(xlimrange )
ax.set_ylim([0, 40] )

plt.xticks(np.arange(180, 210+1, 10.0))
plt.yticks(np.arange(0, 30.0+0.1, 10.0))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_ylim([0, 1] )

fig.savefig('IF1stLast'+'.svg',bbox_inches='tight')