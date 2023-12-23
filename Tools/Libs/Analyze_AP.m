    function [Slope_rise, Slope_decay, Half_width, V_reset_vec, V_thr_vec, V_max_vec]=Analyze_AP(AP_V_mat,dt,O_ind)

% The function is used to analyze the property of an AP matrix.
% The AP_v_mat=[sweep_ind, time_index]. O_index is the index for the max value.
% Slop_rise, Slope_decay, Half_width=[sweep_ind,1];


    Flag_thr=20;  % TO detect the firing threshold.
    
    [n_sweep,n_tind]=size(AP_V_mat);

    AP_dV_mat=(-AP_V_mat(:,1:end-1)+AP_V_mat(:,2:end))/dt;

    Slope_rise= max(AP_dV_mat,[],2);
    Slope_decay= min(AP_dV_mat,[],2);

    % Not sure how I should do the half width. May
    V_max_vec=max(AP_V_mat,[],2);
%     V_min_vec=min(AP_V_mat(:,1:O_ind),[],2);
    
    V_reset_vec= min(AP_V_mat(:,O_ind:end),[],2);
    V_thr_vec=nan(n_sweep,1);
    
    
    V_flag_vec=nan(n_sweep,1);
    
    
    
    Ind_vec=1:n_tind;
    
    Half_width=nan(n_sweep,1);
%     for i=494:n_sweep
    for i=1:n_sweep
        AP_iter=AP_V_mat(i,:);
        AP_dV_iter=(-AP_iter(1:end-1)+AP_iter(2:end))/dt;


%         i
        % For debugging
%         figure(1234)
%         clf
%         hold all
%         plot(Ind_vec*dt,AP_iter)
%         plot(Ind_vec*dt,Ind_vec*0+V_flag_vec(i))
%         plot(Ind1*dt+[0,0],[V_min_vec(i),V_max_vec(i)]  )
%         plot(Ind1*dt+dt1+[0,0],[V_min_vec(i),V_max_vec(i)]  )
%         plot(Ind2*dt+[0,0],[V_min_vec(i),V_max_vec(i)]  )
%         plot(Ind2*dt+dt2+[0,0],[V_min_vec(i),V_max_vec(i)]  )
        
%         i
% 
% 
%         Half_width(i)=(Ind2(1)-Ind1(1))*dt+dt2-dt1;


%         Now adding the detection of firing threshold and the reseting. 
        Flag3=( AP_dV_iter(1:end-1)<=Flag_thr ) & ( AP_dV_iter(2:end)>Flag_thr)  ;
        temp_ind=find(Flag3);
        if isempty(temp_ind)
            continue
        end
        % Using interpolation to get the threshold.
        a1=AP_iter(temp_ind(1)); b1=AP_iter(temp_ind(1)+1);
        a2=AP_dV_iter(temp_ind(1)); b2=AP_dV_iter(temp_ind(1)+1);

        V_thr_vec(i)=b1*(a2-Flag_thr)/(a2-b2)+a1*(b2-Flag_thr)/(b2-a2);
%         V_thr_vec(i)=(AP_iter(temp_ind(1))+AP_iter(temp_ind(1)+1) )/2;

        V_flag_vec(i)=1/2*(V_max_vec(i)+V_thr_vec(i));        
        
        Flag1=( AP_iter(1:O_ind-1)<=V_flag_vec(i) ) & ( AP_iter(2:O_ind)>V_flag_vec(i) )  ;
        Ind1=Ind_vec(Flag1);
        dt1=dt*(V_flag_vec(i)-AP_iter(Ind1)) /(AP_iter(Ind1+1)-AP_iter(Ind1))  ;  
        Flag2=( AP_iter(O_ind:end-1)>=V_flag_vec(i) ) & ( AP_iter(O_ind+1:end)<V_flag_vec(i) )  ;
        Ind2=Ind_vec(Flag2)+O_ind-1;
        dt2=dt*(V_flag_vec(i)-AP_iter(Ind2)) /(AP_iter(Ind2+1)-AP_iter(Ind2));
        if isempty(Ind2) || isempty(Ind1)
            continue
        end
        Half_width(i)=(Ind2(1)-Ind1(1))*dt+dt2-dt1;
    end

end