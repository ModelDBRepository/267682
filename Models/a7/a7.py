# -*- coding: utf-8 -*-
"""
Modified on Sep. 13th, 2021
To mimic the data of a7 17105023
@author: Hongyu Meng
R=136.8*MOhm Timescale=7.22 ms 5 to 15 ms

c=T/R=53*b2.pfarad range from 37 to 110


"""

import brian2 as b2
from datetime import date,datetime
import os
import numpy as np

now=datetime.now()
today=date.today()
current_time=now.strftime("%H%M")
current_day=today.strftime("%m%d_")
timestamp=current_day+current_time

TempPath='SavedResearch/'+timestamp
if os.path.isdir(TempPath)==False:
    os.makedirs(TempPath)
    
    
b2.start_scope()

# Cm=80*pF  gl=5 *nS for I   
delta_T=0.5* b2.mV         # 3 iru 3.5 This one is the original value. I feel I should fit this parameter if I have data.
v_rheobase= -39.0 * b2.mV    # -48 for the sub osc. -60 for the original Bru05 In BerGer05, it's -50.

Inj_max=300*b2.pA
Inj_min=100*b2.pA

C_m=80*b2.pfarad    # Varies a lot. Currently, I choose 50. To decrease the time constant, I decrease the C_m
R_L=120 * b2.Mohm
#g_exp =0 * b2.nS  
#R=500 * b2.Mohm
v_rest=-70.0 * b2.mV    # -60 mV in Litwin code
v_reset=-50.0 * b2.mV
#b= 0 * b2.pA              # 30 for VIP, 80 for the original. 0 for no adaptation.
v_spike=-30. * b2.mV

g_AHP=0.0*b2.nS       # 0.1 for genric model
dw=1            # 30 for VIP, 80 for the original. 0 for no adaptation.



tau_w=500.0 * b2.ms          # default 150 ms
simulation_time=5000 * b2.ms
FIRING_THRESHOLD_v_spike=-10*b2.mV

# Parameters for T-type calcium
gT=8* b2.nS   # 1/5* g_L based on Rinzel 2000
vh= -60 * b2.mV ;v_T= 120 *b2.mV;
tau_hm=5*b2.ms;  # Original 20, reducing that.

tau_hp=100*b2.ms;
      

Num_neu= 200

Inj_vec= b2.linspace(1,Num_neu,Num_neu)/Num_neu *(Inj_max-Inj_min)+Inj_min

# Adding ks current like in StiSej13 or other type
dm= 0.
g_ks=0*b2.nS  # 50* g_L
v_k=-90*b2.mV
tau_s= 200*b2.ms
#group.I='(i+1)/Num_step*Inj_max'                        # Run n simulation to check the f-I curve.
#ds/dt=( 1/(1+exp(0.2*(-53.4-v/b2.mV)  )  ) -  s  )/ tau_s   :1
#ds/dt=( 1/(1+ v/mV )-  s  )/ tau_s    :1



# Oscillating input
A=0.*b2.nA     # Set A=0 to shut down the oscilation input.
f = 10* b2.Hz   # Of course this need to vary, but later.
sigma =10* b2.pA
tau_sigma=0.5*b2.ms

eqs = """
dv/dt = (-(v-v_rest)/R_L +I_exp+I_ks  +  I_stim +I_T +I_adap+I_sigma)/(C_m) : volt (unless refractory)
I_exp=delta_T*exp((v-v_rheobase)/delta_T)/R_L             : amp
dI_sigma/dt = -I_sigma/tau_sigma+ sigma*xi*sqrt(2/tau_sigma)  : amp

I_adap= -g_AHP*w*(v-v_k)       : amp
dw/dt=-w/tau_w : 1

I_ks=-g_ks*mks*(1*h1)*(v-v_k)    :amp
mks_inf=1/(1+exp(-(v/mV+34)/6.5))   :1      # should be 34
dmks/dt= -(mks-mks_inf)/ (6*ms)     :1
hks_inf= 1/(1+exp((v/mV+65)/6.6))   :1
tau_h1= (200+220/(1+exp(-(v/mV+71.6)/6.85)))*ms   : second
dh1/dt= -(h1-hks_inf)/ (tau_h1)     :1

I_stim: amp
I_T= -gT*minf*h*(v- v_T)   :amp
dh/dt= -h/tau_hm *int(v > vh) + (1-h) / tau_hp * int(v< vh)  :1
minf=int(v>vh)  : 1
"""

neuron = b2.NeuronGroup(Num_neu, model=eqs, threshold='v>FIRING_THRESHOLD_v_spike', reset="v=v_reset;w+=dw;mks+=dm", refractory=1.5*b2.ms, method="euler")

# initial values of v and w is set here:
neuron.v = v_rest
neuron.w = 0.0
#neuron.s= 0
#neuron.q= 0
neuron.h1= 0.05

#neuron.I_stim='(i+1)/Num_neu*2*Inj_max-Inj_max'                        # Run n simulation to check the f-I curve.

# Monitoring membrane voltage (v) and w
#statemon= b2.StateMonitor(neuron, ["v", "w","mks","h1","h2","h","minf"], record=True)
statemon= b2.StateMonitor(neuron, ["v", "w","mks","h1"], record=True)

spikemon = b2.SpikeMonitor(neuron, variables='v')

# running simulation
neuron.I_stim=0*b2.pA                       # Run n simulation to check the f-I curve.

b2.run(1000*b2.ms)

#neuron.I_stim='(i+1.0)/Num_neu*Inj_max'                        # Run n simulation to check the f-I curve.
neuron.I_stim='(i+1.0)/Num_neu*(Inj_max-Inj_min)+Inj_min'                        # Run n simulation to check the f-I curve.


b2.run(simulation_time)




#checkregion=b2.array([1000, simulation_time/b2.ms+1000])
#checkregion=b2.array([3000, 4000])

#checkregion=b2.array([0, simulation_time/b2.ms])
#checkregion=b2.array([900, 2900])
checkregion=b2.array([500, 3000])

CheckInd=171  # 19 for 1/2(I_min+I_max)

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
b2.plot(statemon.t/b2.ms, statemon.mks[CheckInd],'C0',label='mks')
b2.plot(statemon.t/b2.ms, statemon.h1[CheckInd],'C2',label='h1')

#b2.plot(statemon.t/b2.ms, statemon.h2[CheckInd],'C1',label='h2')

b2.xlabel('Time (ms)')
b2.ylabel('q (pA)');
#b2.plot(statemon.t/b2.ms, statemon.m[0],'C0',label='m')
#plot(statemon.t/ms, statemon.h2[0],'C4',label='h2')
axes = b2.gca()
axes.set_xlim(checkregion)

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


Adap_ind= b2.zeros((Num_neu,1))
Length_sub= b2.zeros((Num_neu,1))
CVISI_vec= b2.zeros((Num_neu,1))
Std_vec=b2.zeros((Num_neu,1))


spike_trains = spikemon.spike_trains()
for i_neu in range(Num_neu):
    temp=statemon.v[i_neu]/b2.mV
    Std_vec[i_neu]=b2.std(temp[30000::])
    spike_train=spike_trains[i_neu]/b2.ms
    if len(spike_train)>2:
        length=len(spike_train)
        ISI_vec=spike_train[length//2:]-spike_train[length//2-1:-1]
        CVISI_vec[i_neu]=b2.std(ISI_vec)/b2.mean(ISI_vec)
        Adap_ind[i_neu]=1-      (spike_train[1]-spike_train[0])/(spike_train[-1]-spike_train[-2])
        Length_sub[i_neu]=max(ISI_vec)
    

check_vec=b2.array(np.arange(150, 191, 20))


b2.figure(figsize=(4, 4))  # Checking for adaptation index.

for i_neu in check_vec:
    spike_train=spike_trains[i_neu]
    ISI_vec=spike_train[1:]-spike_train[:-1]
    b2.plot((spike_train[1:]+spike_train[:-1])/2,1000/(ISI_vec/b2.ms))

axes = b2.gca()
axes.set_xlim([1, 2])
b2.xlabel('Time (s)')
b2.ylabel('freq (Hz)');
#b2.title(['adaptation index=%.4f' % (1-ISI_vec[0]/ISI_vec[-1])])
b2.savefig(TempPath+'/'+'Freq'+'.jpg')

b2.figure(figsize=(12, 12))  # Checking for adaptation index.

b2.subplot(411)
b2.plot(Inj_vec/b2.nA, spikemon.count/(simulation_time))
b2.xlabel('Injection (nA)')
b2.ylabel('Firing rate (sp/s)')


b2.subplot(412)
b2.plot(Inj_vec/b2.nA,Adap_ind)
b2.xlabel('Injection (nA)')
b2.ylabel('Adap Ind.');

b2.subplot(413)  # plotting the length of subthreshold epoch.
b2.plot(Inj_vec/b2.nA,Std_vec)
b2.xlabel('Injection (nA)')
b2.ylabel('std of 1s to 4s');
ax = b2.gca()

ax.set_ylim([0, 1])


b2.subplot(414)  # plotting the length of subthreshold epoch.
b2.plot(Inj_vec/b2.nA,CVISI_vec)
b2.xlabel('Injection (nA)')
b2.ylabel('CVISI (1)');
b2.savefig(TempPath+'/'+'Summary'+'.jpg')



b2.figure(figsize=(6,6))  # This one is used to generate the phase diagram, currently just one point.

Ind_start=1;
Ind_end=50000;   # 10000 for 1s

volt_trace=statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV
Dif_trace= (statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV-statemon.v[CheckInd,Ind_start-1:Ind_end-1]/b2.mV)/(b2.defaultclock.dt/b2.ms)

b2.plot(volt_trace,Dif_trace)

axes = b2.gca()
axes.set_xlim([-70,-25])
axes.set_ylim([-2, 20])
b2.savefig(TempPath+'/'+('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.nA))+'.jpg')


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
