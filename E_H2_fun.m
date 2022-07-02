function H2 = E_H2_fun(T_air,I_EL,N_EL,n_I_EL,I_FC,N_FC,n_I_FC,E_H2_max,HL_initial)
% H2 = E_H2_fun(T_air,I_EL,N_EL,n_I_EL,I_FC,N_FC,n_I_FC,E_H2_max,HL_initial)
% T is air temperature in C degree.
% I_FC is current of FC in A.
% N_EL is number of EL cells used.
% I_EL is current of EL in A.
% N_FC is number of FC used.
% E_H2_max is H2 tank capacity in kWh.
% E_H2 will be the hydrogrn energy in the tank (kWh).
% HL will be hydrogen level of tank (%).
% P_EL will be energy consumed by electolyzer in kW
% P_FC will be energy produced by fuelcell in kW

% Initiallization
n_inv = 0.95;
E_h2 = HL_initial/100*E_H2_max;      % initial charge of H2 tank in kWh

% Fuel cell and electrolyzer simulation (their power)
[tmp,tmp2,P_EL] = V_EL(I_EL,T_air,N_EL);
[tmp,tmp2,P_FC] = V_FC(I_FC,T_air,N_FC);

% change in total energy of hydrogen system from system perspective
E_h2 = E_h2 + P_EL + P_FC/n_I_FC/n_inv;

if E_h2>E_H2_max
    P_EL = P_EL - (E_h2-E_H2_max);
    E_h2 = E_H2_max;
end

if E_h2<0
    P_FC = P_FC - n_I_FC*n_inv*(E_h2-0);
    E_h2 = 0;
end

% Results
H2.E_H2 = E_h2;
H2.HL = 100*E_h2/E_H2_max;  
H2.E_EL = P_EL/n_I_EL/n_inv;
H2.E_FC = P_FC;

if imag(E_h2)~=0
    error('imag number')
end

end