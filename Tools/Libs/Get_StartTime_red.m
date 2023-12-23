% 08242021 Detect the starting of the Injection 
% I_data: Current data. Should be the maximum sweep.
% Ind_t: Index of the injection started. 
% This one is for red cells, thus we don't have I but V. So changed a bit.
% Comment. It can be used to extract the sag step by using -I_data.
function Ind_t=Get_StartTime_red(I_data)
    I_max=max(I_data);
    I_mid=median(I_data);
    I_flag=(I_max-I_mid)*0.5+I_mid;
    Ind_inj=find(I_data>I_flag);
    Ind_t=Ind_inj(1);

end