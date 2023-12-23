function [Inj_rheo,slope]  = Analyze_SpikeCount(Inj_vec,n_spike_vec,plot_flag, varargin)
    % This code is a simplified version of Analyze_fI that just focus on
    % the ReLU fitting of spike count. 
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
    slope=f.b;
    

    
end