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

Inj_max=300*b2.pA
Inj_min=0*b2.pA

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
simulation_time=5000 * b2.ms
FIRING_THRESHOLD_v_spike=-20*b2.mV

# Parameters for T-type calcium
gT=10* b2.nS   # 2 for the genric model 1/5* g_L based on Rinzel 2000
vh= -60 * b2.mV ;v_T= 120 *b2.mV;
tau_hm=5*b2.ms;   # 20 ms from Rinzel.

tau_hp=100*b2.ms;
      

Num_neu= 300

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
sigma =10* b2.pA
tau_sigma=0.5*b2.ms

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

#neuron.I_stim='(i+1)/Num_neu*2*Inj_max-Inj_max'                        # Run n simulation to check the f-I curve.

# Monitoring membrane voltage (v) and w
#statemon= b2.StateMonitor(neuron, ["v", "w","mks","h1","h2","h","minf"], record=True)
statemon= b2.StateMonitor(neuron, ["v", "w","a","b","I_adap"], record=True)

spikemon = b2.SpikeMonitor(neuron, variables='v')

# running simulation
neuron.I_stim=0                        # Run n simulation to check the f-I curve.

b2.run(1000*b2.ms)

#neuron.I_stim='(i+1.0)/Num_neu*Inj_max'                        # Run n simulation to check the f-I curve.
neuron.I_stim='(i+1.0)/Num_neu*(Inj_max-Inj_min)+Inj_min'                        # Run n simulation to check the f-I curve.


b2.run(simulation_time)

spike_trains = spikemon.spike_trains()

############## plot a7 like initial bursting.

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

checkregion=b2.array([-100, 3000])
#checkregion=b2.array([-10, 100])
CheckInd=120
#CheckInd=236  # 19 for 1/2(I_min+I_max)
#checkregion=b2.array([-100, 1000])



# For better plotting, I need to plot the voltage to the -20 mV

#fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
fig, ax = plt.subplots(dpi=250,figsize=(3,1.7))
spike_train=spike_trains[CheckInd]
Ind_spike=(spike_train/b2.defaultclock.dt).astype(int)
vcopy=statemon.v[CheckInd]/b2.mV
vcopy[Ind_spike]=-20

ax.plot(statemon.t/b2.ms-1000, vcopy, '-k',linewidth=1)   # The injection starts at 0
ax.set_xlim(checkregion)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([-70, -20])

plt.show() 
fig.savefig(('VoltLong%2.0f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')
#fig.savefig(('Volt%2.0f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')




Ind_start=1;
Ind_end=np.int(simulation_time/(b2.defaultclock.dt));   # 10000 for 1s

volt_trace=statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV
Dif_trace= (statemon.v[CheckInd,Ind_start:Ind_end]/b2.mV-statemon.v[CheckInd,Ind_start-1:Ind_end-1]/b2.mV)/(b2.defaultclock.dt/b2.ms)
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))

ax.plot(volt_trace,Dif_trace, '-k',linewidth=1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

ax.set_xlim([-60,-30])
ax.set_ylim([-2, 20])
fig.savefig(('PhasePlanInje%.2f' % (Inj_vec[CheckInd]/b2.pA))+'.svg',bbox_inches='tight')

########## Plot the fI curve  
fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
b2.plot(Inj_vec/b2.pA, spikemon.count/(simulation_time),color='black')
ax.set_xlim([0,300] )
ax.set_ylim([0, 100] )
plt.xticks(np.arange(0, 300+1, 100.0))


#b2.xlabel('Injection (nA)')
b2.ylabel('Total Firing rate (sp/s)')
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_ylim([0, 1] )

fig.savefig('Ficurve'+'.svg',bbox_inches='tight')


spike_trains = spikemon.spike_trains()


fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))


ColorIndex=0
for i_ind in range(6):
    CheckInd=i_ind*40+91
    
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

fig.savefig('1stDifRateTrace'+'.svg',bbox_inches='tight')








#checkregion=b2.array([1000, simulation_time/b2.ms+1000])
#checkregion=b2.array([3000, 4000])

#checkregion=b2.array([0, simulation_time/b2.ms])s
#checkregion=b2.array([900, 2900])
CheckInd=150
checkregion=b2.array([400, 3000])

b2.figure(figsize=(6, 9))
b2.subplot(311)

b2.plot(statemon.t/b2.ms, statemon.v[CheckInd]/b2.mV, '-b')

b2.xlabel('Time (ms)')
b2.ylabel('v (mV)');
b2.title( 'Inj %.4f' % (Inj_vec[CheckInd]/b2.nA) )

axes = b2.gca()
#axes.set_xlim([1000, simulation_time/b2.ms+1000])
axes.set_xlim(checkregion)

#axes.set_xlim([1000, 1000+1000])

b2.subplot(312)
b2.plot(statemon.t/b2.ms, statemon.a[CheckInd],'C0',label='a')
b2.plot(statemon.t/b2.ms, statemon.b[CheckInd],'C2',label='b')

#b2.plot(statemon.t/b2.ms, statemon.h2[CheckInd],'C1',label='h2')

b2.xlabel('Time (ms)')
b2.ylabel('a&b (pA)');
#b2.plot(statemon.t/b2.ms, statemon.m[0],'C0',label='m')
#plot(statemon.t/ms, statemon.h2[0],'C4',label='h2')
axes = b2.gca()
axes.set_xlim(checkregion)
axes.set_ylim([0, 1])

b2.subplot(313)
b2.plot(statemon.t/b2.ms, statemon.I_adap[CheckInd]/b2.pA,'C0',label='a')

#b2.plot(statemon.t/b2.ms, statemon.h2[CheckInd],'C1',label='h2')

b2.xlabel('Time (ms)')
b2.ylabel('I_adap (pA)');
#b2.plot(statemon.t/b2.ms, statemon.m[0],'C0',label='m')
#plot(statemon.t/ms, statemon.h2[0],'C4',label='h2')
axes = b2.gca()
axes.set_xlim(checkregion)


b2.legend()
b2.savefig(TempPath+'/'+('ExampleInj%.4f' % (Inj_vec[CheckInd]/b2.nA))+'.jpg')


##### Now it's time for matching irregular analysis. Also Adaptation


xlimrange=np.array([100 ,300])

CVISI_vec= b2.zeros((Num_neu,1))
Std_vec=b2.zeros((Num_neu,1))

MaxMinratio_vec= b2.zeros((Num_neu,1))
Initrate_vec= b2.zeros((Num_neu,1))
Lastrate_vec= b2.zeros((Num_neu,1))

spike_trains = spikemon.spike_trains()
for i_neu in np.arange(0, Num_neu, 1):
    temp=statemon.v[i_neu]/b2.mV
    Std_vec[i_neu]=b2.std(temp[20000::])
    spike_train=spike_trains[i_neu]/b2.ms
    if len(spike_train)>4:
        length=len(spike_train)
        ISI_vec=spike_train[length//2:-1]-spike_train[(length//2-1):-2]   # Get the latter spike train.
        CVISI_vec[i_neu]=b2.std(ISI_vec)/b2.mean(ISI_vec)
        Initrate_vec[i_neu]=1000/(spike_train[1]-spike_train[0])
        Lastrate_vec[i_neu]=1000/ISI_vec[-1]
        
        ISIratio_temp= np.concatenate([ISI_vec[1:]/ISI_vec[:-1],ISI_vec[:-1]/ISI_vec[1:]])
        MaxMinratio_vec[i_neu]= max(ISIratio_temp)    
        
        
#spike_train=spike_trains[CheckInd]
#ISI_vec=spike_train[1:]-spike_train[:-1]


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
plt.xticks(np.arange(100, 300+1, 100.0))
plt.yticks(np.arange(0, 1.0+0.1, 0.5))

fig.savefig('CVISIMaxmin'+'.svg',bbox_inches='tight')


############# Plot std over injection current  # Maybe same dirty tricks required too.
xlimrange=np.array([100, 300])

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,Std_vec,color='black')
ax.set_xlim([30,150] )

plt.xticks(np.arange(50, 150+1, 50.0))

plt.yticks(np.arange(0, 1.0+0.1, 0.5))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([0, 1] )

fig.savefig('VoltStd'+'.svg',bbox_inches='tight')


############### Now let's do adaptation from PlotFig4

def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def fit_exp_nonlinear(t, y,p0=[-5, 0.01, 25]):
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, t, y, p0, maxfev=10000)
    A, K, C = opt_parms
    return A, K, C




######  Plot different fI curves

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ColorIndex=0
for i_ind in range(6):
    CheckInd=i_ind*40+91
    
    spike_train_all=spike_trains[CheckInd]/b2.ms
    spike_train=spike_train_all[(spike_train_all>100+1000) & (spike_train_all<100+1900) ]   # After initial bursting.

    if  spike_train.size>=2:
        # Non-linear Fit
        ISI_vec=spike_train[1:]-spike_train[:-1]
        FI_vec=1000/(ISI_vec)
        ax.plot(spike_train[:-1]-1000, FI_vec, '.',markersize=3,color=MyColors[ColorIndex,:] )   # The injection starts at 0

        
    if  spike_train.size>=6:
        A, K, C = fit_exp_nonlinear(spike_train[:-1]-1000, 1000/(ISI_vec), [2*FI_vec[-1], -0.001, FI_vec[-1]]) 
        fit_y = model_func(spike_train[:-1]-1000, A, K, C)
        ax.plot(spike_train[:-1]-1000, fit_y, '-',linewidth=1,color=MyColors[ColorIndex,:]) 
    ColorIndex+=1    
#ax.set_xlim([0, 1000])
ax.set_xlim([0, 1000])

ax.set_ylim([0, 100])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(np.arange(0, 1000+1, 500.0))
plt.yticks(np.arange(0, 100+1, 50.0))

fig.savefig('DifRateTrace'+'.svg',bbox_inches='tight')

FitInitRate=np.full(Num_neu,np.nan)     
FitLastRate=np.full(Num_neu,np.nan)     


 
for i_ind in range(Num_neu):
    
    spike_train_all=spike_trains[i_ind]/b2.ms
    spike_train=spike_train_all[(spike_train_all>100+1000) & (spike_train_all<100+2000) ]
#    ISI_vec=spike_train[1:]-spike_train[:-1]

    if  spike_train.size>=20:
        ISI_vec=spike_train[1:]-spike_train[:-1]
        FI_vec=1000/(ISI_vec)
#        ax.plot(spike_train[:-1]-500, 1000/(ISI_vec), '.k',markersize=1)   # The injection starts at 0
        # Non-linear Fit
        A, K, C = fit_exp_nonlinear(spike_train[:-1]-1000, 1000/(ISI_vec), [2*FI_vec[-1], -0.001, FI_vec[-1]])   # last is p0. 
        fit_y = model_func(spike_train[:-1]-1000, A, K, C)
#        ax.plot(spike_train[:-1]-500, fit_y, '-',linewidth=1) 
        FitInitRate[i_ind] =fit_y[1]
        FitLastRate[i_ind] =fit_y[-1]
#    ax.plot(spike_train[:-1]-500, 1000/(ISI_vec), linewidth=1)   # The injection starts at 0



############# Plot 1st IF and last IF over injection current

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,FitInitRate,color='black')
ax.plot(Inj_vec/b2.pA,FitLastRate,color='red')

ax.set_xlim(xlimrange )
#ax.set_ylim([5, 25] )

plt.xticks(np.arange(0, 300+1, 150.0))
#plt.yticks(np.arange(0, 30.0+0.1, 10.0))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_ylim([0, 1] )

fig.savefig('IF1stLast'+'.svg',bbox_inches='tight')


############# Plot adaptation curve

fig, ax = plt.subplots(dpi=250,figsize=(2.2,1.7))
ax.plot(Inj_vec/b2.pA,1-FitLastRate/FitInitRate,color='black')
#ax.plot(Inj_vec/b2.pA,FitLastRate,color='red')

ax.set_xlim(xlimrange )
#ax.set_ylim([5, 25] )

plt.xticks(np.arange(0, 300+1, 150.0))
plt.yticks(np.arange(0, 1+0.01, 0.5))


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([0, 1] )

fig.savefig('AdaptationIndex'+'.svg',bbox_inches='tight')






############## Original analysis
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


b2.savefig(TempPath+'/'+'IrregularAnalyisis'+'.jpg')



from shutil import copyfile      # Save a copy of the code to the saving folder.
copyfile('VIP.py',TempPath+'/'+'VIP.py')

#axes = b2.gca()
#axes.set_xlim([0.0,0.06])

import scipy as sp
from scipy import optimize

def ReLU_func(x,A,B):
    return (x-A)*((x-A)>0).astype(int) *B

def fit_relu(x,y,p0=[100,1]):
    opt_parms, parm_cov = sp.optimize.curve_fit(ReLU_func, x, y, p0, maxfev=10000)
    A, B = opt_parms
    return A, B

Rheo, k=  fit_relu(Inj_vec/b2.pA,spikemon.count/(simulation_time/b2.second),[100,1])   
