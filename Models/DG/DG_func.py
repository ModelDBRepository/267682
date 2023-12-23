# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:50:39 2022
Define the DG function such that can be called iteratively to generate the shuffling gain.
Using the theoretical method to calc the FFT(IAC)
@author: yawnt
"""

import numpy as np
import matplotlib.pyplot as plt
import brian2 as b2
import random
def Shuffle_spikes(spike_train,n):
    # Input: A spike train generated from simulation, n: number of trains that was shuffled.
    # Output: A shuffled spike train.
    CVs= spike_train[1:]-spike_train[:-1]
    CVs= np.insert(CVs,0,spike_train[0])
    output_trains= np.full((n,len(spike_train)),0)
    temp_train=np.full(len(spike_train)+1,0)
    for i in range(n):
        CV_shuffled=CVs
        random.shuffle(CV_shuffled)
        for j in range(len(spike_train)):
            temp_train[j+1]=temp_train[j]+CV_shuffled[j]
        output_trains[i,:]=temp_train[1:]
    return output_trains    
    

def Calc_DG(spike_train_orig,Isig, para, plot_flag=True):
    # Unzip the parameters. All are dimensionless
    timestep= para.timestep
    sim_time= para.simulation_time   # dimensionless
    sigma_I= para.sigma   # dimensionless? Check carefully with the code.
    tau_I= para.tau_sigma     # unit as ms.

    m_window=np.arange(-500.0,500.0,timestep)
    spike_train=spike_train_orig
    spike_train=spike_train[spike_train>1000] # assume first 500ms is the transient period.
    spike_train=spike_train[spike_train<500+sim_time-500] # removing the spiking for the last 500 ms
    
    curr_mat= np.full((len(spike_train),len(m_window)),np.nan) 
    for i_AP in range(len(spike_train)): # for each action potential, truncate a 1s long signal.
        cur_window=spike_train[i_AP]+m_window
        cur_window_ind= (cur_window/timestep).astype(int)
        curr_temp=Isig[cur_window_ind]
        curr_mat[i_AP,:]=curr_temp
        
    curr_ave=np.mean(curr_mat,axis=0)   # autocorrelation of the signal and rate. c_sr  
    
    SPEC_vec=np.arange(0.0,3.001,0.03)  # make sure include the ending point
    spec_vec=10**SPEC_vec
    C_sr_vec=np.full(len(SPEC_vec),np.nan,dtype='complex_')

    for i_freq in range(len(spec_vec)):
        spec_ind= spec_vec[i_freq]
        kernel_vec=   np.exp(-(m_window/1000*spec_ind)**2/2)    # a kernel depends on tau over tau
        cont_sr_vec = curr_ave*kernel_vec*np.exp(-np.pi*spec_ind*m_window/1000*2j)*timestep/1000   # Contribution to C_sr
        C_sr_vec[i_freq]=sum(cont_sr_vec)
 
    C_ss_thry_vec=2*tau_I/1000*sigma_I**2/(1+(2*np.pi*tau_I/1000*spec_vec)**2)
    gain_thry_vec=abs(C_sr_vec)/abs(C_ss_thry_vec)
    phase_thry_vec= np.arctan2(np.imag(C_sr_vec),np.real(C_sr_vec))
    
    
    if plot_flag:
        fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   
        ax.plot(spec_vec, gain_thry_vec,'r',label='thry',linewidth=1)
        ax.set_xlim([1,1000])
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.show()
    
    return gain_thry_vec, phase_thry_vec
    
def Calc_cutoffFreq(gain_thry_vec,plot_flag=True):
    # Calculate the cut_off frequency. Currently, it is 70% of the maximum
    # Also, fit the value from 1 to 80%. 0.8 to 0.64
    SPEC_vec=np.arange(0.0,3.001,0.03)  # make sure include the ending point
    spec_vec=10**SPEC_vec

    gain_norm_vec=gain_thry_vec/gain_thry_vec[0]
    g_vec=np.log(gain_thry_vec/gain_thry_vec[0])/np.log(10)
    ind_vec=np.arange(len(spec_vec)-1,dtype=int)
    
    cut_ind_boolean=(gain_norm_vec[:-1]>=0.7)&(gain_norm_vec[1:]<0.7)
    cut_ind=ind_vec[cut_ind_boolean]
    cut_ind=cut_ind[-1]      # In care there are more than 2 num.
    cut_freq=spec_vec[cut_ind]
    
    fit_ind_boolean=(gain_norm_vec[:-1]>=0.14)&(gain_norm_vec[1:]<0.14)
    fit_ind=ind_vec[fit_ind_boolean]
    fit_ind=fit_ind[0]
    fit_freq=spec_vec[fit_ind]   # In care there are more than 2 num.
    
    gain_fit_vec= g_vec[cut_ind:(fit_ind+1)]
    SPEC_fit_vec= SPEC_vec[cut_ind:(fit_ind+1)]
    spec_fit_vec= spec_vec[cut_ind:(fit_ind+1)]
    ind_fit_vec=ind_vec[cut_ind:(fit_ind+1)]
    
    # maybe I should directly use a log-log fit. Alright, at the weekend.
    
    model=np.polyfit(SPEC_fit_vec,gain_fit_vec,1)   # fit in the log scale.
    k = model[0]
    p= np.poly1d(model)
    fitted_GAIN = p(SPEC_fit_vec)   # In the log space
    fitted_gain = 10**fitted_GAIN
    
    fitted_GAIN_vec= p(SPEC_vec)
    fitted_gain_vec= 10**fitted_GAIN_vec
    # for debug, plot the fitting result.
    if plot_flag:
        fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   
        ax.plot(spec_vec, gain_thry_vec/gain_thry_vec[0],'k',label='exp',linewidth=1)
        ax.plot(spec_fit_vec, fitted_gain,'r',label='fitted',linewidth=1)
        ax.set_xlim([1,1000])
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.show()    
        
        fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   
        ax.plot(SPEC_vec, fitted_GAIN_vec,'k',label='exp',linewidth=1)
        #ax.plot(spec_fit_vec, fitted_gain,'r',label='fitted',linewidth=1)

        plt.show()          
        
        
    return cut_freq, k, fitted_gain_vec   # Cut off frequency, and the fitted slope.
    
