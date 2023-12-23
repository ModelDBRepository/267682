% Calculate the 4th f_I slope and firing rate based on the last 200 ms interval
% I don't believe it will work but Bernardo asked for it, so I need to do it.

function [slope, firing_rate,f_vec]  = FIforth(spike_time_states,Inj_vec,n_spike_vec,Inj_rheo,Inj_slope,plot_flag)

    n_sweeps=length(Inj_vec);

    f_vec=zeros(1,n_sweeps);

    for i=1: n_sweeps
        spike_time=spike_time_states{i};
        spike_time_0= spike_time(spike_time>=800);
        
        if length(spike_time_0)>=2
            ISI_vec_0=spike_time_0(2:end)-spike_time_0(1:end-1);
            f_vec(i)=1000./mean(ISI_vec_0);
        end
        

    end


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
    
    
    % Initialize Fit_ReLU_gain
    Fit_ReLU_gain(0,0,Inj_rheo);
    
    
    ft = fittype( 'Fit_ReLU_gain( x, b )' );
    f = fit( Inj_vec', f_vec', ft, 'StartPoint', Inj_slope,'Exclude', exclude1);

    if plot_flag
        figure(1000);
        clf
        hold all
        plot(f,Inj_vec',f_vec',exclude1)
    end
    
    slope=f.b;
    x=1.7*Inj_rheo;
    
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

    y0=f_vec(Ind_high-1); y1=f_vec(Ind_high);
    firing_rate= y0+(y1-y0)/(x1-x0)*(x-x0);

    if plot_flag
        % varargin(1) should be the filetag
        figure(2000);
        clf
        hold all
        plot(Inj_vec',f_vec')
        plot(x,firing_rate,'o')
     
        
    end

end
