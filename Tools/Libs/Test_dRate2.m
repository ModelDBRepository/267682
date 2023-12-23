function [dRate_score,dRate_vec,dISI_score, dISI_vec, nAP, Score_vec]  = Test_dRate2(spike_time_states,Inj_vec, plot_flag )
% This is the test for ramping behavior

% John: 09132021

    n_sweeps=length(Inj_vec);
    dRate_vec=nan(n_sweeps,1);
    dISI_vec=nan(n_sweeps,1);
    clst_drate=100; clst_ind=[];
    nAP=nan;
    for i=1: n_sweeps
        spike_time=spike_time_states{i};
        Spike_ind_1= spike_time>100  ;  % Get all the spike
        n_spike_1=sum(Spike_ind_1);


        if abs(n_spike_1-9)<clst_drate && n_spike_1>2
            clst_drate=abs(n_spike_1-9);
            clst_ind=i;
            nAP=n_spike_1;
        end        
    
        if sum(Spike_ind_1)>=6 && sum(Spike_ind_1)<=12
            spike_time_stable=spike_time(Spike_ind_1);
            ISI_vec_1=spike_time_stable(2:end)-spike_time_stable(1:end-1);

            x_fit_temp=spike_time_stable(1:n_spike_1-1)';
            x_fit=x_fit_temp-x_fit_temp(1);
%             x_fit=x_fit_temp;

            y_fit_temp=(1000./ISI_vec_1).';
            y_fit=y_fit_temp;
            
            X_fit=[ones(length(x_fit),1) x_fit];    % The shape may be a problem.
            
            b=X_fit\y_fit;
            
%             Ramp_vec(i)=b(2);   % Units are Hz/s
            
            if plot_flag
                figure(i)
                hold on
                plot(x_fit, y_fit,'o-')
                plot(x_fit, x_fit*b(2)+b(1),'r.-')
                hold off
            end
            dRate_vec(i)=b(2)/b(1)*1000;
            dISI_vec(i)= (y_fit(end)-y_fit(1))/y_fit(1);


        end
        
    end

    
    if ~isempty(clst_ind)
        dRate_score=dRate_vec(clst_ind);
        dISI_score=dISI_vec(clst_ind);

        spike_time=spike_time_states{clst_ind};
        Spike_ind_1= spike_time>200  ;  % Get IR score
        Spike_time_stable=spike_time(Spike_ind_1);
        ISI_vec_stable=Spike_time_stable(2:end)-Spike_time_stable(1:end-1);
        CV_ISI= std(ISI_vec_stable)/mean(ISI_vec_stable);
        CV_ISI_score=CV_ISI;  

        rate_ratio_raw=ISI_vec_stable(1:end-1)./ISI_vec_stable(2:end);
        rate_ratio= rate_ratio_raw.*(rate_ratio_raw<=1)+1./rate_ratio_raw.*(rate_ratio_raw>1);
        ISIratio_sweep=1-(rate_ratio);
        ISIratio_score=max(ISIratio_sweep);
        
        if isempty(ISIratio_score)
            ISIratio_score=nan;
        end

        Spike_ind_1= spike_time>100 ;  % Get spike after 100 ms
        n_spike_1=sum(Spike_ind_1);
        if n_spike_1>1
            spike_time_stable=spike_time(Spike_ind_1);
            ISI_vec_1=spike_time_stable(2:end)-spike_time_stable(1:end-1);
    
            x_fit_temp=spike_time_stable(1:n_spike_1-1)';
            x_fit=x_fit_temp-x_fit_temp(1);
    %             x_fit=x_fit_temp;
    
            y_fit_temp=(1000./ISI_vec_1).';
            y_fit=y_fit_temp;
            
            X_fit=[ones(length(x_fit),1) x_fit];    % The shape may be a problem.
            
            b=X_fit\y_fit;

            if plot_flag
                figure(clst_ind)
                hold on
                plot(x_fit, y_fit,'o-')
                plot(x_fit, x_fit*b(2)+b(1),'r.-')
                hold off
            end

            dRate2_score=b(2)*900;
            dISI2_score= (y_fit(end)-y_fit(1))/y_fit(1);
            dISI3_score= y_fit(end)-y_fit(1);
            dISI4_score= dISI3_score/(n_spike_1/0.9);
        else

            dRate2_score=nan;
            dISI2_score= nan;
            dISI3_score=nan;
            dISI4_score=nan;
        end
        
        dInitft_score=1000/(spike_time(2)-spike_time(1));


    else
        dRate_score=nan;
        dISI_score=nan;
        CV_ISI_score=nan;
        ISIratio_score=nan;
        dRate2_score=nan;
        dISI2_score=nan;
        dInitft_score=nan;
        dISI3_score=nan;
        dISI4_score=nan;

    end
    Score_vec=[CV_ISI_score,ISIratio_score,dRate2_score, dISI2_score,dISI3_score,dISI4_score,dInitft_score];
end
