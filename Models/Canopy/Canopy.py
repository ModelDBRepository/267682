# -*- coding: utf-8 -*-
"""
Modified on 06/03/2021
Using the Kv1 from Sciamanna Wilson 2011
To mimic the data of Canopy 16o17034
@author: Hongyu Meng
R=84*MOhm Timescale= 6.5ms
c=T/R=77*b2.pfarad


"""

import brian2 as b2
from datetime import date,datetime
import os
import numpy as np
import matplotlib.pyplot as plt

now=datetime.now()
today=date.today()
current_time=now.strftime("%H%M")
current_day=today.strftime("%m%d_")
timestamp=current_day+current_time

TempPath='SavedResearch/'+timestamp

if os.path.isdir(TempPath)==False:
    os.makedirs(TempPath)


b2.defaultclock.dt = 0.1*b2.ms

b2.start_scope()

# Cm=80*pF  gl=5 *nS for I   
delta_T=1* b2.mV         # 3 iru 3.5 This one is the original value. I feel I should fit this parameter if I have data.
v_rheobase= -38.0 * b2.mV    # -48 for the sub osc. -60 for the original Bru05 In BerGer05, it's -50.

Inj_max=400*b2.pA 
Inj_min=300*b2.pA
C_m=80*b2.pfarad    # Default 77. Increase to match the f-I curve Varies a lot. Currently, I choose 50. To decrease the time constant, I decrease the C_m
R_L=86 * b2.Mohm
#g_exp =0 * b2.nS  
#R=500 * b2.Mohm
v_rest=-67.0 * b2.mV    # -60 mV in Litwin code
v_reset=-52.0 * b2.mV
dw= 0 * b2.pA              # 30 for VIP, 80 for the original. 0 for no adaptation.
v_spike=-30. * b2.mV
tau_w=500.0 * b2.ms          # default 150 ms
simulation_time=10000 * b2.ms
FIRING_THRESHOLD_v_spike=-15*b2.mV

# Parameters for T-type calcium
gT=0* b2.nS   # 1/5* g_L based on Rinzel 2000
vh= -60 * b2.mV ;v_T= 120 *b2.mV;
tau_hm=20*b2.ms; tau_hp=100*b2.ms;
      

Num_neu= 100

Inj_vec= b2.linspace(1,Num_neu,Num_neu)/Num_neu *(Inj_max-Inj_min)+Inj_min 

# Adding ks current like in StiSej13 or other type
dm= 0.

v_shift=20*b2.mV

g_k1= 120*b2.nS
mid_a=-50*b2.mV+v_shift
mid_b=-65*b2.mV+v_shift

sig_a=10*b2.mV
sig_b=-6*b2.mV

tau_a=4*b2.ms      
tau_b=150*b2.ms      # I need to think about this to determine the frequency.

v_k=-90*b2.mV

#group.I='(i+1)/Num_step*Inj_max'                        # Run n simulation to check the f-I curve.
#ds/dt=( 1/(1+exp(0.2*(-53.4-v/b2.mV)  )  ) -  s  )/ tau_s   :1
#ds/dt=( 1/(1+ v/mV )-  s  )/ tau_s    :1



# Oscillating input
sigma =20.* b2.pA

tau_sigma=0.5*b2.ms

eqs = """
dv/dt = (-(v-v_rest)/R_L +I_exp+I_k1  +  I_stim +I_T -  w+I_sigma)/(C_m) : volt (unless refractory)
I_exp=delta_T*exp((v-v_rheobase)/delta_T)/R_L             : amp
dw/dt=-w/tau_w : amp

dI_sigma/dt = -I_sigma/tau_sigma+ sigma*xi*sqrt(2/tau_sigma)  : amp

I_k1=-g_k1*(a*a*a*b)*(v-v_k)  : amp
da/dt = (ainf-a)/tau_a  :   1
db/dt = (binf-b)/tau_b  :   1         (unless refractory)
ainf=1/(1+exp(-(v- mid_a) /sig_a      )   )   :1
binf=1/(1+exp(-(v- mid_b) /sig_b      )   )   :1

I_stim: amp
I_T= -gT*minf*h*(v- v_T)   :amp
dh/dt= -h/tau_hm *int(v > vh) + (1-h) / tau_hp * int(v< vh)  :1
minf=int(v>vh)  : 1
"""


neuron = b2.NeuronGroup(Num_neu, model=eqs, threshold='v>FIRING_THRESHOLD_v_spike', reset="v=v_reset;w+=dw;b+=-0.002", refractory=2*b2.ms, method="euler")

# initial values of v and w is set here:
neuron.v = v_rest
neuron.w = 0.0 * b2.pA
#neuron.s= 0
#neuron.q= 0
neuron.b= 0.3

#neuron.I_stim='(i+1)/Num_neu*2*Inj_max-Inj_max'                        # Run n simulation to check the f-I curve.

# Monitoring membrane voltage (v) and w
#statemon= b2.StateMonitor(neuron, ["v", "w","mks","h1","h2","h","minf"], record=True)
statemon= b2.StateMonitor(neuron, ["v", "w","a","b"], record=True)

spikemon = b2.SpikeMonitor(neuron, variables='v')

# running simulation
neuron.I_stim=0                        # Run n simulation to check the f-I curve.

b2.run(500*b2.ms)

#neuron.I_stim='(i+1.0)/Num_neu*Inj_max'                        # Run n simulation to check the f-I curve.
neuron.I_stim='(i+1.0)/Num_neu*(Inj_max-Inj_min)+Inj_min'                        # Run n simulation to check the f-I curve.


b2.run(simulation_time)




#checkregion=b2.array([1000, simulation_time/b2.ms+1000])
#checkregion=b2.array([3000, 4000])

#checkregion=b2.array([0, simulation_time/b2.ms])
#checkregion=b2.array([900, 2900])
checkregion=b2.array([0, 3000])
CheckInd=79  # 19 for 1/2(I_min+I_max)

b2.figure(figsize=(12, 12))
b2.subplot(211)

b2.plot(statemon.t/b2.ms, statemon.v[CheckInd]/b2.mV, '-b')

b2.xlabel('Time (ms)')
b2.ylabel('v (mV)');
b2.title( 'Inj %.4f' % (Inj_vec[CheckInd]/b2.nA) )

axes = b2.gca()
#axes.set_xlim([1000, simulation_time/b2.ms+1000])
axes.set_xlim(checkregion)

#axes.set_xlim([1000, 1000+1000])

b2.subplot(212)
b2.plot(statemon.t/b2.ms, statemon.a[CheckInd],'C0',label='a')
b2.plot(statemon.t/b2.ms, statemon.b[CheckInd],'C2',label='b')

#b2.plot(statemon.t/b2.ms, statemon.h2[CheckInd],'C1',label='h2')

b2.xlabel('Time (ms)')
b2.ylabel('a & b');
#b2.plot(statemon.t/b2.ms, statemon.m[0],'C0',label='m')
#plot(statemon.t/ms, statemon.h2[0],'C4',label='h2')
axes = b2.gca()
axes.set_xlim(checkregion)
axes.set_ylim([0.2, 0.3])

b2.legend()
b2.savefig(TempPath+'/'+('ExampleInj%.4f' % (Inj_vec[CheckInd]/b2.nA))+'.svg')

#b2.subplot(313)
#b2.plot(spikemon.t/b2.ms, spikemon.v[0]/b2.mV, 'ob')
#b2.legend()

#b2.xlabel('s ')
#b2.ylabel('v (mV)')



#b2.plot(spikemon.t/b2.ms, spikemon.v/b2.mV, 'ob')

#db2.plot( statemon.s[0],statemon.v[0]/b2.mV,'C3',label='s')
#plot( statemon.n[0],statemon.v[0]/mV,'C1',label='n')
#b2.plot( statemon.ht[0],statemon.v[0]/b2.mV,'C2',label='ht')

'''
b2.subplot(414)
b2.plot(statemon.t/b2.ms, statemon.h[CheckInd],'C0',label='h')
b2.plot(statemon.t/b2.ms, statemon.minf[CheckInd],'C2',label='minf')
b2.legend()
b2.xlabel('Time (ms)')
b2.ylabel('h, minf (1)')
#b2.plot(statemon.t/b2.ms, statemon.m[0],'C0',label='m')
#plot(statemon.t/ms, statemon.h2[0],'C4',label='h2')
axes = b2.gca()
axes.set_xlim(checkregion)
'''
#avda_vec=np.full(v_num,np.nan)


Adap_ind= b2.zeros((Num_neu,1))
Length_sub= b2.zeros((Num_neu,1))
CVISI_vec= np.full((Num_neu,1),np.nan)  # b2.zeros((Num_neu,1))
Std_vec=b2.zeros((Num_neu,1))

MaxMinratio_vec= b2.zeros((Num_neu,1))
Initrate_vec= b2.zeros((Num_neu,1))
Lastrate_vec= b2.zeros((Num_neu,1))

spike_trains = spikemon.spike_trains()
for i_neu in range(Num_neu):
    temp=statemon.v[i_neu]/b2.mV
    Std_vec[i_neu]=b2.std(temp[10000::])
    spike_train=spike_trains[i_neu]/b2.ms
    if len(spike_train)>5:
        length=len(spike_train)
        ISI_vec=spike_train[length//2:-1]-spike_train[(length//2-1):-2]
        CVISI_vec[i_neu]=b2.std(ISI_vec)/b2.mean(ISI_vec)
        Adap_ind[i_neu]=1-      (spike_train[1]-spike_train[0])/(spike_train[-1]-spike_train[-2])
        Length_sub[i_neu]=max(ISI_vec)
        Initrate_vec[i_neu]=1000/(spike_train[1]-spike_train[0])
        Lastrate_vec[i_neu]=1000/ISI_vec[-1]
        
        ISIratio_temp= np.concatenate([ISI_vec[1:]/ISI_vec[:-1],ISI_vec[:-1]/ISI_vec[1:]])
        MaxMinratio_vec[i_neu]= max(ISIratio_temp)    
        
        
spike_train=spike_trains[CheckInd]
ISI_vec=spike_train[1:]-spike_train[:-1]




b2.figure(figsize=(12, 8))  # Checking for adaptation index.

b2.subplot(212)



b2.plot((spike_train[1:]+spike_train[:-1])/2,1000/(ISI_vec/b2.ms))
b2.xlabel('Time (ms)')
b2.ylabel('freq (Hz)');
#b2.title(['adaptation index=%.4f' % (1-ISI_vec[0]/ISI_vec[-1])])

b2.subplot(211)
b2.hist(ISI_vec/b2.ms,50)
b2.xlabel('ISI (ms)')
b2.ylabel('Counts ');
b2.title( 'CV ISI= %.4f' % (b2.std(ISI_vec)/b2.mean(ISI_vec)) )


b2.savefig(TempPath+'/'+'CVISI'+'.svg')

b2.figure(figsize=(12, 12))  # Checking for adaptation index.

b2.subplot(411)
b2.plot(Inj_vec/b2.nA, spikemon.count/(simulation_time))
b2.xlabel('Injection (nA)')
b2.ylabel('Firing rate (sp/s)')
b2.title('fig 2')

b2.subplot(412)
b2.plot(Inj_vec/b2.nA,Adap_ind)
b2.xlabel('Injection (nA)')
b2.ylabel('Adap Ind.');

b2.subplot(413)  # plotting the length of subthreshold epoch.
b2.plot(Inj_vec/b2.nA,Std_vec)
b2.xlabel('Injection (nA)')
b2.ylabel('std of 1s to 4s');
axes = b2.gca()
axes.set_ylim([0., 1])

b2.subplot(414)  # plotting the length of subthreshold epoch.
b2.plot(Inj_vec/b2.nA,CVISI_vec)
b2.xlabel('Injection (nA)')
b2.ylabel('CVISI (1)');
b2.savefig(TempPath+'/'+'Summary'+'.svg')

b2.figure(figsize=(12, 12))  # Checking for adaptation index.

b2.subplot(411)
b2.plot(Inj_vec/b2.nA, spikemon.count/(simulation_time))
#b2.xlabel('Injection (nA)')
b2.ylabel('Firing rate (sp/s)')


b2.subplot(412)
b2.plot(Inj_vec/b2.nA,MaxMinratio_vec)
#b2.xlabel('Injection (nA)')
b2.ylabel('Max/Min ISI (1)');

b2.subplot(413)  # plotting the length of subthreshold epoch.
b2.plot(Inj_vec/b2.nA,CVISI_vec)
#b2.xlabel('Injection (nA)')
b2.ylabel('CVISI (1)');

b2.subplot(414)  # plotting the length of subthreshold epoch.
b2.plot(Inj_vec/b2.nA,Initrate_vec)
b2.plot(Inj_vec/b2.nA,Lastrate_vec)
b2.xlabel('Injection (nA)')
b2.ylabel('Init & Last  (Hz)');


b2.savefig(TempPath+'/'+'IrregularAnalyisis'+'.svg')

b2.figure(figsize=(6,6))  # This one is used to generate the phase diagram, currently just one point.

Ind_start=1;
Ind_end=np.int(simulation_time/(b2.defaultclock.dt));   # 10000 for 1s

volt_trace=statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV
Dif_trace= (statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV-statemon.v[CheckInd,Ind_start-1:Ind_end-1]/b2.mV)/(b2.defaultclock.dt/b2.ms)

b2.plot(volt_trace,Dif_trace)

axes = b2.gca()
axes.set_xlim([-55,-25])
axes.set_ylim([-8, 13])
b2.savefig(TempPath+'/'+('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.nA))+'.svg')


'''
b2.figure(figsize=(12, 12))  # Checking for adaptation index.
b2.subplot(211)
b2.plot(statemon.mks[CheckInd], statemon.v[CheckInd]/b2.mV,'C0')

b2.ylabel('v (mV)')
b2.xlabel('mks (1)');
b2.subplot(212)
b2.plot( statemon.h1[CheckInd], statemon.v[CheckInd]/b2.mV,'C0')

b2.ylabel('v (mV)')
b2.xlabel('h1 (1)');
'''
from shutil import copyfile      # Save a copy of the code to the saving folder.
copyfile('Canopy.py',TempPath+'/'+'Canopy.py')

#
#axes = b2.gca()
#axes.set_xlim([0.0,0.06])
#%% Try to get the ReLu fitting from the model.
import scipy as sp
from scipy import optimize

def ReLU_func(x,A,B):
    return (x-A)*((x-A)>0).astype(int) *B

def fit_relu(x,y,p0=[100,1]):
    opt_parms, parm_cov = sp.optimize.curve_fit(ReLU_func, x, y, p0, maxfev=10000)
    A, B = opt_parms
    return A, B

Rheo, k=  fit_relu(Inj_vec/b2.pA,spikemon.count/(simulation_time/b2.second),[100,1])   

#b2.plot(Inj_vec/b2.nA, spikemon.count/(simulation_time))



print('\a')

