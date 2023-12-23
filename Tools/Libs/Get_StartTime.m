% 08242021 Detect the starting of the Injection 
% I_data: Current data. Should be the maximum sweep.
% Ind_t: Index of the injection started. 

% Comment. It can be used to extract the sag step by using -I_data.
function [Ind_t, Ind_tend]=Get_StartTime(I_data)
    I_max=max(I_data);
    I_mid=median(I_data);
    I_flag=(I_max-I_mid)*0.8+I_mid;
    Ind_inj=find(I_data>I_flag);
    Ind_t=Ind_inj(1);
    Ind_tend=Ind_inj(end);
end