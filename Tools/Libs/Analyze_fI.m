function [Inj_rheo,slopes, firing_rates,f_vec1,f_vec2, f_vec3, figlab1]  = Analyze_fI(spike_time_states,Inj_vec,n_spike_vec,plot_flag, varargin)
    % This code is used to analyze firing rate that general 4 rates.
    % Spike count, average initial firing rate, fitted rate 
    % Output: rheobase. Defined by the spike count
    % fI_vec1, 2, 3: Spike count, average initial firing rate, fitted rate over inj vec, average last 200 ms
    % fI_1, 2, 3: Measured firing rate at 1.8 rheobase. Calculated by interpolation or extrapolation
    % Slope_1, 2, 3: Slopes measured by reLR firring, Maybe just fit gain.
    
    % Notice, firing rate 2 may be negative because of fitting and long delay of the 1st spike
    %---------
    
    
    ft = fittype( 'Fit_ReLU( x, a, b )' );
    n_sweeps=length(Inj_vec);
    fit_flag=ones(1,n_sweeps);
    fir_ind=find(n_spike_vec>0);
    if length(fir_ind)>5 
        fit_flag(1:fir_ind(5))=0;
    else
        fit_flag=0;
    end
    
    [~,max_fir_ind]=max(n_spike_vec);
    
    exc_fir_ind=(1:n_sweeps) > max_fir_ind;
    
    exclude1 = (n_spike_vec> 40) & fit_flag | exc_fir_ind  ;   % Fitting around the rheobase
    f = fit( Inj_vec', n_spike_vec', ft, 'StartPoint', [min(Inj_vec ),0.5],'Exclude', exclude1);
    
    if plot_flag
        figlab1=figure(1000);
        clf
        hold all
        plot(f,Inj_vec',n_spike_vec',exclude1)
    end
    
    Inj_rheo=f.a;
    
    %  Generate next two fI curve
    n_sweeps=length(Inj_vec);
    f_vec1=zeros(1,n_sweeps);
    f_vec2=zeros(1,n_sweeps);
    f_vec3=zeros(1,n_sweeps);

    for i=1: n_sweeps
        spike_time=spike_time_states{i};
        spike_time1= spike_time(spike_time<200);
        spike_time2= spike_time(spike_time<500);
        spike_time3= spike_time(spike_time>=800);
        if length(spike_time1)>=2
            ISI_vec1=spike_time1(2:end)-spike_time1(1:end-1);
            f_vec1(i)=1000./mean(ISI_vec1);
        end
        
        if length(spike_time2)>=4
            ISI_vec2=spike_time2(2:end)-spike_time2(1:end-1);
            InsRate_vec=1000./ISI_vec2;
            p=polyfit(spike_time2(1:end-1), InsRate_vec , 2);
            f_vec2(i)=polyval(p,0);
            
            % For debugging
%             figure(1011)
%             hold on
%             plot(spike_time2(1:end-1),InsRate_vec,'o')
%             f_test=polyval(p,spike_time2(1:end-1));
%             plot(spike_time2(1:end-1),f_test,'r--')
%             
%             i
        end        
        
        if length(spike_time3)>=2
            ISI_vec3=spike_time3(2:end)-spike_time3(1:end-1);
            f_vec3(i)=1000./mean(ISI_vec3);
        end
                
        
    end
    

    % Initialize Fit_ReLU_gain
    Fit_ReLU_gain(0,0,Inj_rheo);
    
    ft = fittype( 'Fit_ReLU_gain( x, b )' );
    f1 = fit( Inj_vec', f_vec1', ft, 'StartPoint', [f.b],'Exclude', exclude1);
    
    f2 = fit( Inj_vec', f_vec2', ft, 'StartPoint', [f.b],'Exclude', exclude1);
    f3 = fit( Inj_vec', f_vec3', ft, 'StartPoint', [f.b],'Exclude', exclude1);
    
    if plot_flag
            % varargin(1) should be the filetag
        fig1=figure(1000);
        plot(f1,Inj_vec',f_vec1',exclude1)
        plot(f2,Inj_vec',f_vec2',exclude1)
        plot(f3,Inj_vec',f_vec3',exclude1)
        figname= 'FittedfIcurve';
        file_name= strcat(varargin{1},'_',figname,'.fig');
        hgsave(fig1,file_name) 
    end    

    slopes=[f.b,f1.b,f2.b,f3.b];
    
    
    
    x=1.7*Inj_rheo;
    
    % Now let's do the interpolation and extrapolation.

    if x<max(Inj_vec)
        Ind_temp= find(x<Inj_vec);
        Ind_high=Ind_temp(1);
    else
        Ind_high=n_sweeps;
    end
    
    if Ind_high==1
        Ind_high=2;
    end
        
    
    x0=Inj_vec(Ind_high-1); x1=Inj_vec(Ind_high);

    y0=n_spike_vec(Ind_high-1); y1=n_spike_vec(Ind_high);
    firing_rate0= y0+(y1-y0)/(x1-x0)*(x-x0);

    y0=f_vec1(Ind_high-1); y1=f_vec1(Ind_high);
    firing_rate1= y0+(y1-y0)/(x1-x0)*(x-x0);

    y0=f_vec2(Ind_high-1); y1=f_vec2(Ind_high);
    firing_rate2= y0+(y1-y0)/(x1-x0)*(x-x0);        
        
    y0=f_vec3(Ind_high-1); y1=f_vec3(Ind_high);
    firing_rate3= y0+(y1-y0)/(x1-x0)*(x-x0); 

    
%     %Alternatively, let's fit a quadratic function, then do the measurement
%     Ind_nonzero=n_spike_vec>0;
%     x_vec= Inj_vec(Ind_nonzero);
%     y0_vec=n_spike_vec(Ind_nonzero);
%     y1_vec=f_vec1(Ind_nonzero);
%     y2_vec=f_vec2(Ind_nonzero);
%     
%     p0=polyfit(x_vec,y0_vec,2);
%     p1=polyfit(x_vec,y1_vec,2);
%     p2=polyfit(x_vec,y2_vec,2);
%     
%     firing_rate0=polyval(p0,x);
%     firing_rate1=polyval(p1,x);
%     firing_rate2=polyval(p2,x);
    
    firing_rates=[firing_rate0,firing_rate1,firing_rate2,firing_rate3];
    if plot_flag
        % varargin(1) should be the filetag
        fig2=figure(2000);
        clf
        hold all
        plot(Inj_vec',n_spike_vec')
        plot(Inj_vec',f_vec1')
        plot(Inj_vec',f_vec2')
        plot(Inj_vec',f_vec3')

        plot([x,x,x,x],firing_rates,'o')
        figname= 'FirngRates';
        file_name= strcat(varargin{1},'_',figname,'.fig');
        hgsave(fig2,file_name)        
        
    end
    
end