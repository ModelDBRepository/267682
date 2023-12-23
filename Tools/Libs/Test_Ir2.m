function [Flag,ISIratio_vec,CV_ISI_vec  ]  = Test_Ir2(spike_time_states,Inj_vec, MaxMin_thr, ISI_thr, plot_flag )
%   A test for the irregularity. Here is based on the consequencial two
%   ISIs.

    if plot_flag
        close all
    end
    n_sweeps=length(Inj_vec);
    
    ISIratio_vec=nan(n_sweeps,1);

    CV_ISI_vec=nan(n_sweeps,1);
    
    for i=1: n_sweeps
        
        spike_time=spike_time_states{i};

        
        Spike_ind_1= spike_time>(500)  ;  % Get the steady part
        Spike_time_stable=spike_time(Spike_ind_1);
        if sum(Spike_ind_1)>=5 && length(spike_time)<=50
            ISI_vec_stable=Spike_time_stable(2:end)-Spike_time_stable(1:end-1);


            CV_ISI= std(ISI_vec_stable)/mean(ISI_vec_stable);
            CV_ISI_vec(i)=CV_ISI;  

% Using the maximum ratio between the subsequent ISI to calculate the ISIratio            
            rate_ratio_raw=ISI_vec_stable(1:end-1)./ISI_vec_stable(2:end);
            rate_ratio= min(rate_ratio_raw,1./rate_ratio_raw);
%             rate_ratio=(max(ISI_vec_stable)-min(ISI_vec_stable))/max(ISI_vec_stable);
            ISIratio_vec(i)=1-min(rate_ratio);
            
            if plot_flag
                figure(100+i)   % For debug
                clf
                plot(Spike_time_stable(1:end-1),ISI_vec_stable,'.k','markersize',15)      
                set(gcf,'units','points','position',[0+i*15,200,400,300])
                title(sprintf('CVISI=%2.3g,MaxminRatio=%2.3g', CV_ISI,rate_ratio ))
            end
            
            
                        
        end


    end
    
    Flag1=nanmax(ISIratio_vec)>MaxMin_thr;
    Flag2=nanmax(CV_ISI_vec)>ISI_thr;
%     Flag= Flag1 | Flag2;
    Flag=  Flag1;
end