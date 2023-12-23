% This script is used to test the data Shlomo generated
close all
addpath('./Libs/');

Flag_sag=0;
Iter_plot_flag=true;
% load('Info_Cells')
% ind_cell=51-1;
% load('L1Info.mat')
CType='a7';
if strcmp(CType,'NGFC') 
    window_t=5;
elseif strcmp(CType,'Canopy')
    window_t=4;
%     timescale_0=1/5
else
    window_t=4;
end
    
file_tag='18503040';     % T-type adaptation


[d,si,h]=abfload(strcat(file_tag,'.abf'));

v_data=squeeze(d(:,1,:));   
I_data=squeeze(d(:,2,:));      
% abfload for p-clamp


[t_len,swipe_length]=size(v_data);
checkind_vec=1:swipe_length;      % The index of swipes

% checkind_vec=11;

dt=0.05;    % ms

% The range of recording. [ms]
ind_max=min(t_len,round(10000/0.05)-1);

ind=1:ind_max;
% ind=round(800/0.05):(round(1200/0.05)-1);
lind=length(ind);
t_axis=dt*(1:t_len);

% Paras for figures
xInit=50;    yInit=350;
width=320;    height=225;


FFTind=round(400/0.05):(round(1300/0.05)-1);    % Works for inj. < rheobase
% OUind=round(500/0.05):(round(1200/0.05)-1);    % Works for inj. < rheobase
OUind=round(400/0.05):(round(1300/0.05)-1);    % Works for inj. < rheobase

% Paras for fft
Fs = 1000/dt;            % Sampling frequency                    
T = dt;             %dT
L = length(FFTind);             % Length of signal, ms
t = (0:L-1)*dt/1000;        % Time vector
% S = 0.7*sin(2*pi*100*t) + sin(2*pi*400*t);   fourier analysis sanity test
% X = S + 0*randn(size(t));
% y= fft(X);

% Allocation for output. n_spike_vec=zeros(1,max(checkind_vec));
Adap_vec=zeros(1,max(checkind_vec));
Adap_vec2=zeros(1,max(checkind_vec));   % Adaptation Index, AI, 1c


AdapTau_vec=zeros(1,max(checkind_vec));  % Timescale of Adaptation, AI, 1c

OU_mean_vec=zeros(1,max(checkind_vec));
OU_sigma_vec=zeros(1,max(checkind_vec));


Inj_vec=zeros(1,max(checkind_vec));
CV_ISI_vec= zeros(1,max(checkind_vec));
rate_ratio_vec= nan(1,max(checkind_vec));   % (Max-min)/Max of ISI,   1d
frate_vec=  zeros(1,max(checkind_vec));      % Steady firing rate, 1a
n_spike_vec=  zeros(1,max(checkind_vec));   

ISI_init_vec= nan(1,max(checkind_vec));
ISI_last_vec= nan(1,max(checkind_vec));

delay_vec=nan(1,max(checkind_vec));

timescale_vec=nan(1,max(checkind_vec));
R_vec=nan(1,max(checkind_vec));
Capacity_vec=nan(1,max(checkind_vec));

t_inj_ind=Get_StartTime(I_data(:,max(checkind_vec)));  
t_inj=t_axis(t_inj_ind);

I_baseline=mean(I_data(  round((t_inj-50)/dt):round(t_inj/dt), max(checkind_vec) ));


[t_sag_ind,t_sagend_ind]=Get_StartTime(-I_data(:,max(checkind_vec)));  % This one is used to get the start/end time of the sag.

t_sag=t_axis(t_sag_ind);
t_sagend=t_axis(t_sagend_ind);

Fir_thr=Get_SpikeThres(v_data(ind,max(checkind_vec)));

v_rest= median(v_data(1:round(t_inj_ind-100/dt),max(checkind_vec)));


shift_t=1;
FstAPmat=[];
ScdAPmat=[];
LatterAPmat=[];
LatterLowSpikeAPmat=[];


states= cell(swipe_length,1);      % A variable that save everything. 
indexes= nan(1, max(checkind_vec));   % A matrix that save all kinds of readouts.  Need a doc to explain details.
    % (1,:)  Injection vector.
    
temp_max_n_spike=0;    
% checkind_vec=9
for checkind=checkind_vec
    
% for checkind=[1:23,25:max(checkind_vec)]
% for checkind=40
    
Sample_trace=v_data(ind,checkind); 

if Iter_plot_flag
    figure(checkind*10+1)
    plot(t_axis(ind)-t_inj,Sample_trace,'k','linewidth',1)


    xlim([-100,1000])
    ylim([-75,-20])
    box off
    set(gca,'Ytick',-70:10:-10)

    set(gcf,'units','points','position',[checkind*8,yInit,width,height])
    figname= sprintf('SampleTraceInj%2.0f',Inj_vec(checkind)-I_baseline);
    ax=gca;
    ax.FontSize=16;
%     file_name= strcat(fileloc,figname);
%     saveas(figlab1,file_name,'emf')   

%     set(gcf,'units','points','position',[xInit+checkind*15,yInit,width,height])
end


% Find the spike and calc the spike number.

flag=(Sample_trace(1:lind-1)<Fir_thr) & (Sample_trace(2:lind)>=Fir_thr);


% Not sure why there are spikes before injection and after. Clear after
% that.

[ind_spike]=find(flag);


% Exlude the data points over 2000 Hz. or 10 data points. 
min_spike_int= 10;
exclude_ind_flag= (ind_spike(2:end)-ind_spike(1:end-1))<10;
ind_spike([false;exclude_ind_flag])=[];


spike_temp=t_axis(ind(ind_spike));
spike_time  =spike_temp( (spike_temp>=t_inj) &  (spike_temp<t_inj+1001));  

states{checkind}=spike_time-t_inj;


n_spike=length(spike_time);   % Maybe can shrink this out. Not for now. 
n_spike_vec(checkind)=n_spike;

spike_ind=round(spike_time/dt);

OU_mean_vec(checkind)= mean(Sample_trace(t_inj_ind+round(500/dt):t_inj_ind+round(900/dt)    ));
OU_sigma_vec(checkind)= std(Sample_trace(t_inj_ind+round(500/dt):t_inj_ind+round(900/dt)    ));

% Need a function to calculate ISI, CV_ISI
if n_spike>0
    delay_vec(checkind)=spike_time(1)-t_inj;

end


if n_spike>1
    ISI_vec=spike_time(2:n_spike)-spike_time(1:n_spike-1);

    % Only using spike_time>200 ms to cal CV_ISI to avoid the impact from
    % adaptation.
%     spike_time_stable=spike_time(spike_time>spike_time(1)+100);   % THis is used to avoid initial bursting
    spike_time_stable=spike_time( (spike_time>t_inj+500) &  (spike_time<t_inj+1001));     % Looking at the firing rate at the latter half. 
    if length(spike_time_stable)>=2
        ISI_vec_stable=spike_time_stable(2:end)-spike_time_stable(1:end-1);
        mean_ISI=mean(ISI_vec_stable);

        frate_vec(checkind)=1000/mean_ISI;

        Std_ISI = std(ISI_vec_stable);
        CV_ISI= Std_ISI/mean_ISI;
        CV_ISI_vec(checkind)=CV_ISI;    
        rate_ratio=(max(ISI_vec_stable)-min(ISI_vec_stable))/max(ISI_vec_stable);
        rate_ratio_vec(checkind)=rate_ratio;
    end


    if Iter_plot_flag
        figlab1= figure(checkind*10+2);  % Plot some phase diagram


        Ind_zoom=round((t_inj-50)/0.05):(round((t_inj+1000)/0.05)-1);
        t_axis_zoom=t_axis(Ind_zoom);
        L_Ind_zoom=length(Ind_zoom);
        slope_mat=nan(L_Ind_zoom,1);
        v_mat=nan(L_Ind_zoom,1);
    
        smooth_window=round(0.25/0.05); % Using 0.25 ms as the smooth resolute
        Ind_smooth=min(Ind_zoom):smooth_window:max(Ind_zoom); % may have a dt diff.
        t_axis_smooth=t_axis(Ind_smooth);
        l_smooth=length(t_axis_smooth);    
        
        temp_trace=v_data(Ind_zoom,checkind);
        temp_mat=reshape(temp_trace,smooth_window,round(L_Ind_zoom/smooth_window));
        Smooth_trace=mean(temp_mat,1);    

        plot(Smooth_trace(2:end-1)',(Smooth_trace(3:end)'-Smooth_trace(1:end-2)')/(smooth_window*dt),'k','linewidth',1)
        xlim([-55,-30])
        ylim([-2,10])
        box off
        set(gcf,'units','points','position',[xInit+checkind*8,yInit-200,width,height])
%         figname= sprintf('DvDtInj%2.0f',Inj_vec(checkind)-I_baseline);
        ax=gca;
        ax.FontSize=20;
    end





    if Iter_plot_flag
        figlab1=figure(checkind*10+3);
    end
    
    if n_spike<=5
        if Iter_plot_flag

            plot(spike_time(1:n_spike-1)-t_inj,1000./ISI_vec,'o')
            title('Freq [Hz] over time')
            set(gcf,'units','points','position',[xInit+checkind*8,yInit-400,width,height]) 
        end
    % Now try to fit a curve. Only fit for the curve if n>20  Seems that if
    % include the first datapoint, there will be a problem. So starting at
    % the second datapoint
    else    % Fitting to generate the adaptation index
        x_fit_temp=spike_time(2:n_spike-1)';
        x_fit=x_fit_temp-t_inj;
        y_fit_temp=(1000./ISI_vec(2:end)).';
        y_fit=y_fit_temp/y_fit_temp(1);
        g_model = fittype('a-b*exp(-c*x)');

        f0_fitted = fit(x_fit,y_fit,g_model,'StartPoint',[0;-10;0.00]);
        xx = linspace(0,1000,100);

        if Iter_plot_flag
            hold all
      
            plot(spike_time(1:n_spike-1)-t_inj,1000./ISI_vec,'ko','MarkerFaceColor','k')
            plot(xx,f0_fitted(xx)*y_fit_temp(1),'k-','linewidth',3);
%             title(['Freq [Hz] over time, Adap Ind=',num2str(-f0_fitted.b),'Last/Init',num2str(f0_fitted(xx(end))/f0_fitted(xx(1))),'Exp fact.=',num2str(f0_fitted.c)]);
            set(gcf,'units','points','position',[checkind*8,yInit-200,width,height]) 
            ax=gca;
            ax.FontSize=16;
            ylim([0, inf])
        end
        
        Adap_vec(checkind)=-f0_fitted.b;
        Adap_vec2(checkind)=1-f0_fitted(xx(end))/f0_fitted(xx(1));
        AdapTau_vec(checkind)=1/f0_fitted.c;

        
        
    end
    ISI_init_vec(checkind)= ISI_vec(1);
    ISI_last_vec(checkind)= ISI_vec(end);
end



Injind=round((t_inj+500)/0.05):(round((t_inj+900)/0.05)-1);  % Automatic pick the time to check the injection current
Inj_I= mean(I_data(Injind,checkind));
Inj_vec(checkind)=Inj_I;     % ?pA


[FstAPiter, ScdAPiter,LatterAPiter] = Get_APs(Sample_trace, spike_ind, shift_t,window_t);

FstAPmat=[FstAPmat;FstAPiter];
ScdAPmat=[ScdAPmat;ScdAPiter];

LatterAPmat=[LatterAPmat;LatterAPiter];
if n_spike<=40 && n_spike>temp_max_n_spike  % The latter one is to exclude the cases where the cell has reached the max firing rate.
    LatterLowSpikeAPmat=[LatterLowSpikeAPmat;LatterAPmat];
end

temp_max_n_spike=max(temp_max_n_spike,n_spike);

fitInd=round((t_sagend+1)/0.05):round((t_sagend+101)/0.05);   % Here is after releasing. 
% fitInd=round(2300/0.05):(round(2400/0.05)-1);   
% % fitInd=round(8510/0.05):(round(8600/0.05)-1);
% % fitInd=round(1750/0.05):(round(1800/0.05)-1);
% % fitInd=round(8400/0.05):(round(8500/0.05)-1);
x_tau= t_axis(fitInd)'-t_axis(fitInd(1));
y_tau=-(v_data(fitInd,checkind)-v_data(fitInd(end),checkind));   
g_model = fittype('a+b*exp(-c*x)');
% For fitting debug

%     figlab2=figure(checkind*10+2);
%     plot(x_tau,y_tau)

f0= fit(x_tau,y_tau,g_model,'StartPoint',[0;10;0.1]);

if 1/f0.c<100 && 1/f0.c>0
    timescale_vec(checkind)=1/f0.c;
end

if Iter_plot_flag && Flag_sag
    figlab2=figure(checkind*10+2);
    plot(f0,x_tau,y_tau)
    title(sprintf('fitting of r, start time= %2.3g', t_axis(fitInd(1)) ))
    
    set(gcf,'units','points','position',[xInit+checkind*15,yInit-200,width,height]) 

end



% % Using the data at 2280 to 2290 to calc the I_, V. 2130 to 2140 to cal the
% % I, V before.

depInd=round((t_sagend-20)/0.05):round((t_sagend-10)/0.05);
nomInd=round((t_sag-20)/0.05):round((t_sag-10)/0.05);


% depInd=round(2280/0.05):(round(2290/0.05)-1);   % Here is during sagging
% 
% % depInd=round(8380/0.05):(round(8390/0.05)-1);
% % nomInd=round(8230/0.05):(round(8240/0.05)-1);
% 
% 
V_dep=mean(v_data(depInd,checkind)); I_dep=mean(I_data(depInd,checkind));
V_nom=mean(v_data(nomInd,checkind)); I_nom=mean(I_data(nomInd,checkind));
R_temp=(V_nom-V_dep)/(I_nom-I_dep)*1000;  % Units: mV/pA*10^3=  Mohm
R_vec(checkind)=R_temp;                     % He
Capacity_vec(checkind)=timescale_vec(checkind)/R_temp;




end

indexes(1,:) = Inj_vec;

% For output
% file_loc=['Saved_Research/',file_tag];

% I_data=squeeze(d(:,Icha,:));      % What's recorded here? Why such large noise?

figure(5000)   % Check the recording of channels 
c_ind=max(checkind);
subplot(2,1,1)
plot(t_axis(ind),v_data(ind,c_ind))
title('Voltage channel')
subplot(2,1,2)
plot(t_axis(ind),I_data(ind,c_ind))

[Inj_rheo,slopes, firing_rates,f_vec1,f_vec2, ~]  = Analyze_fI(states,Inj_vec,n_spike_vec,1,file_tag);

figlab1=figure(1000);   % with figlab1




figname= 'FiringRates';
file_name= strcat(file_tag,'_',figname,'.fig');
hgsave(figlab1,file_name)        

% zoom_range=1:7;  



figlab4=figure(4000);
clf
subplot(2,3,1)
plot(Inj_vec,n_spike_vec,'linewidth',3)
title('Num Spike over Inj current')

subplot(2,3,2)
hold all
plot(Inj_vec,Adap_vec,'linewidth',3)
plot(Inj_vec,Adap_vec2,'linewidth',3)

title('Adaptation Index over Inj current')
ylim([-0.1 2])
subplot(2,3,3)
plot(Inj_vec,AdapTau_vec,'linewidth',3)
title('Adaptation Tau[ms] over Inj current')
ylim([-10 1000])

subplot(2,3,4)
plot(Inj_vec,OU_mean_vec,'linewidth',3)
title('mean [mV], 0.4s to 1.3s')
subplot(2,3,5)
plot(Inj_vec,OU_sigma_vec,'linewidth',3)
title('std [mV]')
subplot(2,3,6)
hold all
plot(Inj_vec,1000./ISI_init_vec,'linewidth',3)
plot(Inj_vec,1000./ISI_last_vec,'linewidth',3)
hold off


set(gcf,'units','points','position',[xInit,yInit-height,width*3,height*1])   % Putting the figure on the left screen. You may not see one with only 1 screen!

title('Init&last firing rate')
figname= 'Summary';
file_name= strcat(file_tag,'_',figname);
saveas(figlab4,file_name,'emf')        
saveas(figlab4,file_name,'fig')        




figlab6=figure(6001);   % This one is for generating the comparable stats for NGFC cells.
clf
subplot(4,1,1)
plot(Inj_vec,n_spike_vec,'linewidth',3)
title('Num Spike over Inj current')
ax=gca;
ax.FontSize=16;
subplot(4,1,2)
plot(Inj_vec,rate_ratio_vec,'.-','linewidth',3)

ax=gca;
ax.FontSize=16;
title('Max/min of firing rate')

% xlim([150 250])

set(gcf,'units','points','position',[xInit,yInit,width*3,height*1])

subplot(4,1,3)
plot(Inj_vec,CV_ISI_vec,'linewidth',3)
ax=gca;
ax.FontSize=16;
title('CV ISI over Inj current')

subplot(4,1,4)
hold all
plot(Inj_vec,1000./ISI_init_vec,'linewidth',3)
plot(Inj_vec,1000./ISI_last_vec,'linewidth',3)
hold off

title('Init&last firing rate [Hz]')
ax=gca;
ax.FontSize=16;
% xlim([150 250])
set(gcf,'units','points','position',[xInit+100,yInit-200,width,height*1.5])
figname= 'Irregular Analysis';
file_name= strcat(file_tag,'_',figname);
saveas(figlab6,file_name,'fig')        
saveas(figlab6,file_name,'emf')   


% Plotting the APs.
TimeVec=-shift_t:dt:window_t-shift_t-dt;

[Slope_rise1, Slope_decay1, Half_width1]=Analyze_AP(FstAPmat,dt,round(shift_t/dt)+1);
[Slope_rise2, Slope_decay2, Half_width2]=Analyze_AP(ScdAPmat,dt,round(shift_t/dt)+1);

figlab8=figure(8000);
clf
subplot(3,2,1)
% ylim([-60 50])
plot(TimeVec,FstAPmat)

title(sprintf('1st AP Distribution, mean HW= %2.3g', nanmean(Half_width1)))

subplot(3,2,3)
if ~isempty(ScdAPmat)
    plot(TimeVec,ScdAPmat)
    title(sprintf('2nd AP Distribution, mean HW= %2.3g', nanmean(Half_width2)))

end
% ylim([-60 50])


subplot(3,2,5)
% ylim([-60 50])
if ~isempty(LatterLowSpikeAPmat)

    plot(TimeVec,LatterLowSpikeAPmat)

    title('Latter Low AP Distribution')
end


subplot(3,2,2)
AP_dV_mat=(-FstAPmat(:,1:end-1)+FstAPmat(:,2:end))/dt;
AP_V_mat=(FstAPmat(:,1:end-1)+FstAPmat(:,2:end))/2;

% ylim([-60 50])
plot(AP_V_mat',AP_dV_mat'       )

title('1st AP dV')

subplot(3,2,4)
if ~isempty(ScdAPmat)

% ylim([-60 50])
AP_dV_mat=(-ScdAPmat(:,1:end-1)+ScdAPmat(:,2:end))/dt;
AP_V_mat=(ScdAPmat(:,1:end-1)+ScdAPmat(:,2:end))/2;

% ylim([-60 50])
plot(AP_V_mat',AP_dV_mat'       )
end

title('2nd AP dV')


subplot(3,2,6)
if ~isempty(LatterLowSpikeAPmat)

AP_dV_mat=(-LatterLowSpikeAPmat(:,1:end-1)+LatterLowSpikeAPmat(:,2:end))/dt;
AP_V_mat=(LatterLowSpikeAPmat(:,1:end-1)+LatterLowSpikeAPmat(:,2:end))/2;

plot(AP_V_mat',AP_dV_mat'        )
% ylim([-60 50])

title('latter AP dV')
end
figname= 'AP analysis';
file_name= strcat(file_tag,'_',figname);
saveas(figlab8,file_name,'fig')        
saveas(figlab8,file_name,'emf')   
set(gcf,'units','points','position',[xInit+100,yInit-400,width,height*2])

% Get the min and max


figure(9000)

% [Slope_rise, Slope_decay, Half_width, V_reset_vec, V_thr_vec, V_max_vec]=Analyze_AP(LatterAPmat,dt,round(shift_t/dt)+1);
[Slope_rise, Slope_decay, Half_width, V_reset_vec, V_thr_vec, V_max_vec]=Analyze_AP(LatterLowSpikeAPmat ,dt,round(shift_t/dt)+1);  % Using this to avoid adapt in the form of AP

subplot(4,1,1)
plot(Slope_rise)
title(sprintf('Latter APs, Ave Rise Slope %2.3g',mean(Slope_rise)))

subplot(4,1,2)

plot(Slope_decay)
title(sprintf('Ave Decay Slope %2.3g',mean(Slope_decay)))

subplot(4,1,3)
plot(Half_width)
title(sprintf('Ave Half Width %2.3g',mean(Half_width)))

subplot(4,1,4)
hold all
plot(V_reset_vec)
plot(V_thr_vec)
plot(V_max_vec)

title(sprintf('Reset=%2.4g, Threshold=%2.4g, Max V=%2.4g',mean(V_reset_vec),mean(V_thr_vec), mean(V_max_vec)  ))





% [Slope_rise, Slope_decay, Half_width,V_reset_vec,V_thr_vec, V_max_vec]=Analyze_AP(LatterAPmat,dt,round(shift_t/dt)+1);

fig10=figure(10000); % Summarize of the new stuff.
clf
subplot(6,1,1)

plot(Inj_vec,delay_vec)

delaynotnan=delay_vec(~isnan(delay_vec));
delay_1st=delaynotnan(1);
title((sprintf('first Delay= %2.3g',delay_1st)))

subplot(6,1,2)
hold all
plot(Inj_vec,n_spike_vec)   % Last 500 ms

% f_ind=find(frate_vec>0);
% 
% % Matlab linear regression on the non-zero datapoints
% X=[ones(length(f_ind),1), Inj_vec(f_ind)'];
% Adap_slope=X\frate_vec(f_ind)';
% 
% % Fit a relu function and plot
% ft = fittype( 'Fit_ReLU( x, a, b )' );
% 
% fit_flag=ones(1,swipe_length);
%     if swipe_length>5 
%         fit_flag(1:5)=0;
%     else
%         fit_flag=0;
%     end
% exclude1 = (n_spike_vec> 40) & fit_flag  ;   % Fitting around the rheobase
% f = fit( Inj_vec', frate_vec', ft, 'StartPoint', [min(Inj_vec )-10,0.5],'Exclude', exclude1);
% plot(f,Inj_vec',frate_vec',exclude1)
% 
% title((sprintf('f-I Slope %2.3g, ReLU fit, rheobase= %2.3g, slope= %2.3g',Adap_slope(2), f.a, f.b)))
% 
% 

subplot(6,1,3)
hold all
plot(Inj_vec,timescale_vec)   % Last 500 ms
title(sprintf('Ave Timescale=%2.3g',mean(timescale_vec(~isnan(timescale_vec))) ))

subplot(6,1,4)
hold all

plot(Inj_vec,R_vec)   % Last 500 ms
title(sprintf('Ave R=%2.3g',mean(R_vec) ))


subplot(6,1,5)

plot(Inj_vec,Capacity_vec)   % Last 500 ms
title(sprintf('Ave C=%2.3g',mean(Capacity_vec) ))



subplot(6,1,6)
plot(Inj_vec,rate_ratio_vec)
title(sprintf('(Max-Min)/Max mean=%2.4g',max(rate_ratio_vec)))





set(gcf,'units','points','position',[xInit+800,yInit*0.1,width,height*2])

figname= '3DAnalysis';
file_name= strcat(file_tag,'_',figname,'.fig');
hgsave(fig10,file_name)        

% value_output=[Adap_vec2(end),Adap_slope(2), f.a, f.b, Slope_rise2(1),Slope_decay2(1),Half_width2(1),mean(V_reset_vec),mean(V_thr_vec), mean(V_max_vec),max(rate_ratio_vec)];

% This Value Output is arranged as: Reset	Threshold	Half Width (latter) Delay, TImescale, REsistance,	Capacity	Rheobase	fISlope    Rise Slope 	Decay Slope

value_output=[mean(V_reset_vec),nanmean(V_thr_vec), nanmean(Half_width),delay_1st,nanmean(timescale_vec),mean(R_vec),nanmean(Capacity_vec)*1000,  Inj_rheo, slopes(1), slopes(2),slopes(3),slopes(4), mean(Slope_rise),   mean(Slope_decay) ];
    

% value_output=[Adap_vec2(end),Half_width2(1), Inj_rheo,firing_rates(1),firing_rates(2),firing_rates(3), slopes(1), slopes(2),slopes(3), Inj_rheo*1.7/max(Inj_vec)];

open('value_output')

% Saving all the stuff here.


    figname= 'struct';
    file_name= strcat(file_tag,'_',figname,'.mat');
    save(file_name, '-regexp', '^(?!(d|si|h|v_data|I_data|fig10|figlab1|figlab3|figlab4|figlab6|figlab8|L1Info)$).');
% 
% figure(11000)  % This is for steady firing rate
% clf
% hold all
% plot(Inj_vec,frate_vec)
% 
% f_ind=find(frate_vec>0);
% 
% % Matlab linear regression on the non-zero datapoints
% X=[ones(length(f_ind),1), Inj_vec(f_ind)'];
% Adap_slope=X\frate_vec(f_ind)';
% 
% % Fit a relu function and plot
% ft = fittype( 'Fit_ReLU( x, a, b )' );
% fit_flag=ones(1,swipe_length);
% fit_flag(1:5)=0;
% 
% exclude1 = (n_spike_vec> 40) & fit_flag  ;   % Fitting around the rheobase
% f = fit( Inj_vec', frate_vec', ft, 'StartPoint', [min(Inj_vec ),0.5],'Exclude', exclude1);
% plot(f,Inj_vec',frate_vec',exclude1)
% 
% title((sprintf('f-I Slope %2.3g, ReLU fit, rheobase= %2.3g, slope= %2.3g',Adap_slope(2), f.a, f.b)))
