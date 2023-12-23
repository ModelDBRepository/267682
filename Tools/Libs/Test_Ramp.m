function [Flag,Ramp_vec, Delay_1st ]  = Test_Ramp(spike_time_states,Inj_vec, Ramp_thr,  plot_flag )
% This is the test for ramping behavior

% John: 09132021

    n_sweeps=length(Inj_vec);
    Ramp_vec=nan(n_sweeps,1);
    Delay_vec=nan(n_sweeps,1);

    for i=1: n_sweeps
        spike_time=spike_time_states{i};
        Spike_ind_1= spike_time>0  ;  % Get all the spike
        
        if sum(Spike_ind_1)>=3
            spike_time_stable=spike_time(Spike_ind_1);
            n_spike_1=sum(Spike_ind_1);
            ISI_vec_1=spike_time_stable(2:end)-spike_time_stable(1:end-1);

            x_fit_temp=spike_time_stable(1:n_spike_1-1)';
            x_fit=x_fit_temp-x_fit_temp(1);
%             x_fit=x_fit_temp;

            y_fit_temp=(1000./ISI_vec_1).';
            y_fit=y_fit_temp;
            
            X_fit=[ones(length(x_fit),1) x_fit];    % The shape may be a problem.
            
            b=X_fit\y_fit;
            
            Ramp_vec(i)=b(2);   % Units are Hz/s
            

            figure(i)
            plot(x_fit, y_fit,'o-')
            
        end
        
        if sum(Spike_ind_1)>=2
             Delay_vec(i)=spike_time(Spike_ind_1(1));
        end

    end


    Flag1=nanmin(Ramp_vec)>Ramp_thr;   % Usually should be 0
    
    Delay_notnan=Delay_vec(~isnan(Delay_vec));
    if ~isempty(Delay_notnan)
         Delay_1st=Delay_notnan(1);
    else
        Delay_1st=nan;
    end
    Flag= Flag1 ;    
end
