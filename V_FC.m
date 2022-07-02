function [Vfc,HCR,P_fc] = V_FC(I_FC,T_FC,N_FC)
% [Vfc,HCR,P_fc] = V_FC(I_FC,T_FC,N_FC)
% I_FC is current of FC in A.
% N_FC is number of FC cells.
% T_FC is considered to be equal to air temperature (as an estimate).
% V_FC will be voltage of FC in V.
% HCR will be hydrogen consumption rate in L/hr (in this project equals to
% L).
% P_el will be power of FC in kW.

% Initillization
T = size(I_FC,2);
Vfc = zeros(1,T);
n_inv = 0.95;
Vbus = 48; % in V

% FC properties
V_FC0 = 33.18;              % FC constant in V
C1_FC = -0.013;             % FC constant in V/C
C2_FC = -1.57;              % FC constant in V
I_FC0 = 8.798;              % FC constant in A
R_FC = -2.04;               % FC constant in ohm.C
C_H2 = 2.39;                % conversion factor in Ah/L
n_I_FC = 0.45;              % FC efficiency
% in Kelouwani model hydrogen pressure is considered to be 1 atm, while in
% fact we use compressures, the energy conversion stayes true.
dH = 286;                   % Conversion constant (kJ/mol)
V_T = 22.4;                 % Conversion constant (L/mol)

% max. current limit of 35  PEM-FC cell
if abs(I_FC)/N_FC>(7e3/n_I_FC/Vbus/35)
    error('Over discharge of FC');
    I_FC = -(7e3/n_I_FC/Vbus)*fix(N_FC/35);
end

% because of the log function in Vfc, I_FC must be positive
I_FC = -I_FC;

% Calculation
Vfc = V_FC0*ones(1,T) + C1_FC*T_FC + C2_FC*log(I_FC/I_FC0) + R_FC*I_FC./T_FC;
math_error_indx = Vfc==inf;
Vfc(math_error_indx) = 0;

HCR = N_FC*(n_I_FC/C_H2)*(I_FC*Vbus/Vfc*n_inv);
HCR(math_error_indx) = 0;
P_fc = -HCR*dH/V_T/3600;          % FC power in kW
end