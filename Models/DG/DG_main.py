# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 17:37:15 2022
Using this us a script that call different DG_func to calculate DG remained
properties. Should run a long simulation first to get the data first.
@author: yawnt
"""

import DG_func as df
import brian2 as b2
import numpy as np
import matplotlib.pyplot as plt
import time

spike_trains = spikemon.spike_trains()
spike_train_orig=spike_train=spike_trains[0]/b2.ms
Isig=statemon.I_sigma[0]/b2.pA

class para:
    pass
para.timestep=b2.defaultclock.dt/b2.ms
para.simulation_time= simulation_time/b2.ms    # dimensionless
para.sigma=sigma/b2.pA     # dimensionless? Check carefully with the code.
para.tau_sigma=tau_sigma/b2.ms     # unit as ms.

gain_thry_vec= df.Calc_DG(spike_train_orig,Isig, para)

n_shuffle=500
shuffled_trains=df.Shuffle_spikes(spike_train_orig,n_shuffle)
gain_shuffle_mat= np.full((n_shuffle,len(gain_thry_vec)),np.nan)
t=time.time()
for i_shuffle in range(n_shuffle):
    if (i_shuffle % 100 ==0): 
        print('Num of shuffle', i_shuffle)
    gain_shuffle_mat[i_shuffle,:]= df.Calc_DG(shuffled_trains[i_shuffle,:],Isig, para, plot_flag=False)
elapsed=time.time()-t

gain_shuffle_ave= np.mean(gain_shuffle_mat,axis=0)
gain_shuffle_std= np.std(gain_shuffle_mat,axis=0)
conf_cutoff= 1.645   # norm.ppf(0.95)
gain_conf_vec=gain_shuffle_ave+ conf_cutoff*gain_shuffle_std

cut_freq, k, fitted_gain_vec = df.Calc_cutoffFreq(gain_thry_vec,plot_flag=False)

# plotting here
SPEC_vec=np.arange(0.0,3.001,0.03)  # make sure include the ending point
spec_vec=10**SPEC_vec

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   
ax.plot(spec_vec, gain_thry_vec/gain_thry_vec[0],'k',label='exp',linewidth=2)
ax.plot(spec_vec, gain_conf_vec/gain_thry_vec[0],'silver',label='shuffle',linewidth=1)
ax.plot(spec_vec[5:], fitted_gain_vec[5:],'--',color='silver',label='fitted',linewidth=1)
ax.set_xlim([1,1000])
ax.set_ylim([0.01,2])
ax.set_xscale('log')
ax.set_yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show() 

cellname='VIP' 
firingrate=spikemon.count[0]*1000*b2.ms/(simulation_time+500*b2.ms)
fig.savefig(cellname+('%.1f'% firingrate)+('cutoff%.3f_slope%.3f_sigma%.1f_Inj%.2f' % (cut_freq,k,para.sigma,Inj_max/b2.pA))+'.svg',bbox_inches='tight')
print('\a')
