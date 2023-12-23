function [FstAPmat, ScdAPmat, LatterAPmat] = Get_APs(V_data, spike_ind,shift_t,window_t)
    % Get APs of a trace. 
    % The default value used for this function. I may change these into persistant variables if it's needed.
    dt=0.05;    % ms
    Ind_wind=round(window_t/dt);  % The window to find a spike.
           % Max time is at shift_t

    % Not sure this works.. Let's see
    n_spike=length(spike_ind);
    
    if n_spike>0 

        FstAPmat=nan(1,Ind_wind);
    else
        FstAPmat=[];
    end
    
    if n_spike>1 

        ScdAPmat=nan(1,Ind_wind);
    else
        ScdAPmat=[];
    end
    if n_spike>2
        LatterAPmat=nan(n_spike-2,Ind_wind);
    else
        LatterAPmat=[];
    end

    for i_spike=1: n_spike
        % search the max volt in a spike wave window.
        Cur_ind=spike_ind(i_spike);

        TempVolt=V_data(Cur_ind:Cur_ind+Ind_wind);
        [~,Ind_max]=max(TempVolt);

        Tuned_ind=Cur_ind+Ind_max-1;
        TempWaveForm=V_data((Tuned_ind-round(shift_t/dt)):(Tuned_ind+round(  (window_t-shift_t)  /dt)-1));
        if i_spike==1

            FstAPmat(:)=TempWaveForm;
        elseif i_spike==2
            ScdAPmat(:)=TempWaveForm;
        else
            LatterAPmat(i_spike-2,:)=TempWaveForm; % Seperate the waveform of the first spike or the latter spikes.
        end
    end



end

