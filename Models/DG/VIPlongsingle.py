# -*- coding: utf-8 -*-
"""
Modified on 06/25/2021
Using the Kv1 from Sciamanna Wilson 2011
Challenge myself to mimic the VIP cell  10052021
@author: Hongyu Meng
R=225*MOhm Timescale= 10.3ms
c=T/R=46.3*b2.pfarad


"""

import brian2 as b2
from datetime import date,datetime
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

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
v_rheobase= -40.0 * b2.mV    # -48 for the sub osc. -60 for the original Bru05 In BerGer05, it's -50.

Inj_max=25*b2.pA
Inj_min=25*b2.pA
sigma =50* b2.pA

C_m=80*b2.pfarad    # Varies a lot. Currently, I choose 50. To decrease the time constant, I decrease the C_m
R_L=400 * b2.Mohm
#g_exp =0 * b2.nS  
#R=500 * b2.Mohm
v_rest=-69.0 * b2.mV    # -60 mV in Litwin code
v_reset=-50.0 * b2.mV

g_AHP=0.1 *b2.nS       # 0.1 for genric model
dw=1            # 30 for VIP, 80 for the original. 0 for no adaptation.
v_spike=-30. * b2.mV
tau_w=500.0 * b2.ms          # default 150 ms
simulation_time=1000000 * b2.ms
FIRING_THRESHOLD_v_spike=-20*b2.mV

# Parameters for T-type calcium
gT=10* b2.nS   # 2 for the genric model 1/5* g_L based on Rinzel 2000
vh= -60 * b2.mV ;v_T= 120 *b2.mV;
tau_hm=5*b2.ms;   # 20 ms from Rinzel.

tau_hp=100*b2.ms;
      

Num_neu= 1

Inj_vec= b2.linspace(1,Num_neu,Num_neu)/Num_neu *(Inj_max-Inj_min)+Inj_min

# Adding ks current like in StiSej13 or other type
dm= 0.

v_shift=20*b2.mV

g_k1= 150*b2.nS      #150 is good to introduce the behavior I want
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
A=0.*b2.nA     # Set A=0 to shut down the oscilation input.
f = 10* b2.Hz   # Of course this need to vary, but later.
tau_sigma=25*b2.ms

eqs = """
dv/dt = (-(v-v_rest)/R_L +I_exp+I_k1  +  I_stim +I_T + I_adap+I_sigma)/(C_m) : volt (unless refractory)
I_exp=delta_T*exp((v-v_rheobase)/delta_T)/R_L             : amp

I_adap= -g_AHP*w*(v-v_k)       : amp
dw/dt=-w/tau_w : 1

dI_sigma/dt = -I_sigma/tau_sigma+ sigma*xi*sqrt(2/tau_sigma)  : amp

I_k1=-g_k1*(a*a*a*b)*(v-v_k)  : amp
da/dt = (ainf-a)/tau_a  :   1
db/dt = (binf-b)/tau_b  :   1    (unless refractory)
ainf=1/(1+exp(-(v- mid_a) /sig_a      )   )   :1
binf=1/(1+exp(-(v- mid_b) /sig_b      )   )   :1

I_stim: amp
I_T= -gT*minf*h*(v- v_T)   :amp
dh/dt= -h/tau_hm *int(v > vh) + (1-h) / tau_hp * int(v< vh)  :1
minf=int(v>vh)  : 1
"""

neuron = b2.NeuronGroup(Num_neu, model=eqs, threshold='v>FIRING_THRESHOLD_v_spike', reset="v=v_reset;w+=dw;b+=-0.001", refractory=2*b2.ms, method="euler")

# initial values of v and w is set here:
neuron.v = v_rest
neuron.w = 0.0 
#neuron.s= 0
#neuron.q= 0
neuron.b= 0.3
neuron.h= 0

#neuron.I_stim='(i+1)/Num_neu*2*Inj_max-Inj_max'                        # Run n simulation to check the f-I curve.

# Monitoring membrane voltage (v) and w
#statemon= b2.StateMonitor(neuron, ["v", "w","mks","h1","h2","h","minf"], record=True)
statemon= b2.StateMonitor(neuron, ["v", "w","a","b","I_sigma"], record=True)

spikemon = b2.SpikeMonitor(neuron, variables='v')

# running simulation
#neuron.I_stim=0                        # Run n simulation to check the f-I curve.

#b2.run(500*b2.ms)

#neuron.I_stim='(i+1.0)/Num_neu*Inj_max'                        # Run n simulation to check the f-I curve.
neuron.I_stim='(i+1.0)/Num_neu*(Inj_max-Inj_min)+Inj_min'                        # Run n simulation to check the f-I curve.


b2.run(simulation_time+500*b2.ms, report='text')



Vstd_vec= np.full(Num_neu,np.nan) 
Istd_vec= np.full(Num_neu,np.nan) 

for i_sweep in range(Num_neu):
    Vstd_vec[i_sweep]=np.std(statemon.v[i_sweep,5000:]/b2.mV)
    Istd_vec[i_sweep]=np.std(statemon.I_sigma[i_sweep,5000:]/b2.pA)

b2.figure(figsize=(12, 12))
b2.subplot(311)
b2.plot(Inj_vec/b2.nA, spikemon.count/(simulation_time+500*b2.ms))
b2.xlabel('Inj (pA)')
b2.ylabel('nAPs');
b2.title( 'fI')

b2.subplot(312)
b2.plot(Inj_vec/b2.nA, Vstd_vec)
b2.xlabel('Inj (pA)')
b2.ylabel('volt std (mV)');
b2.subplot(313)
b2.plot(Inj_vec/b2.nA, Istd_vec)
b2.xlabel('Inj (pA)')
b2.ylabel('Isigma std (pA)');



checkregion=b2.array([5000, 5500])
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))   
ax.plot(statemon.t/b2.ms, statemon.v[0]/b2.mV, '-k',linewidth=0.2)
#ax.xlabel('Time (ms)')
#ax.ylabel('v (mV)');
ax.set_xlim(checkregion)


print('\a')

#b2.plot(Inj_vec/b2.nA, spikemon.count/(simulation_time))

