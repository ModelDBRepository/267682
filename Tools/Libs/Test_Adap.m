% This is the test for initial bursting.
% John: 09132021

function [Flag,Adap_vec2,Timescale_vec  ]  = Test_Adap(spike_time_states,Inj_vec, Adap_thr, Timescale_thr, plot_flag )
% Input form: spike_time_states,Inj_vec, Adap_thr, Timescale_thr, plot_flag

    n_sweeps=length(Inj_vec);
    Adap_vec=nan(n_sweeps,1);
    Adap_vec2=nan(n_sweeps,1);

    Timescale_vec=nan(n_sweeps,1);

    for i=1: n_sweeps
        spike_time=spike_time_states{i};
        Spike_ind_1= spike_time>(100)  ;  % Get the latter spiking part
        
        if sum(Spike_ind_1)>=6
            spike_time_1=spike_time(Spike_ind_1);
            n_spike_1=sum(Spike_ind_1);
            ISI_vec_1=spike_time_1(2:end)-spike_time_1(1:end-1);

            x_fit_temp=spike_time_1(1:n_spike_1-1)';
            x_fit=x_fit_temp-x_fit_temp(1);
%             x_fit=x_fit_temp;

            y_fit_temp=(1000./ISI_vec_1).';
            y_fit=y_fit_temp/y_fit_temp(1);
            g_model = fittype('a+b*exp(-c*x)');
%             
%             figure(1)   % For debug
%             clf
%             hold all
%             plot(x_fit,y_fit,'ok')
            
            f0_fitted = fit(x_fit,y_fit,g_model,'StartPoint',[0;1;0.]);
            xx = linspace(min(x_fit),max(x_fit),100);

            Adap_vec(i)=-f0_fitted.b;
            Adap_vec2(i)=1-f0_fitted(xx(end))/f0_fitted(xx(1));
            Timescale_vec(i)=1/f0_fitted.c;
            
            if plot_flag
                 plot(f0_fitted,x_fit',y_fit')
            end
        end
        
        

    end


%     Flag1=nanmin(Adap_vec2)>Adap_thr;
%     Flag2=nanmax(Timescale_vec)<Timescale_thr;
%     Flag= Flag1 && Flag2;

    Flag1=nanmin(Adap_vec2)>Adap_thr;
%     Flag2=nanmax(Timescale_vec)<Timescale_thr;
    Flag= Flag1 ;    
end
