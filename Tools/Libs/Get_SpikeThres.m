function V_spike_thr = Get_SpikeThres(V_data)
%Generate an automated threshold to DETECT spikes. 
    V_max=max(V_data);
    V_max=min(V_max,0);  % Added here because some a7 can go as high to 40
    V_mid=median(V_data);
    V_spike_thr=max(V_max*0.8+V_mid*0.2,-20);  % Changed from 0.6 to 0.8 for cell 44
    
        
    
end

