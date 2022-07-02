function [Vel,HGR,P_el] = V_EL(I_EL,T_EL,N_EL)
% [Vel,HGR,P_el] = V_EL(I_EL,T_EL,N_EL)
% I_EL is current of EL in A.
% N_EL is number of EL cells.
% T_EL is considered to be equal to air temperature (as an estimate).
% V_EL will be voltage of EL in V.
% HGR is hydrogen generation rate in L/hr (in this project equals to L).
% P_el will be power of EL in kW.

% Initiallization
T = size(I_EL,2);
Vbus = 48; % in V
n_inv = 0.95;
% Vel = zeros(1,T);

% EL properties
V_EL0 = 22.25;              % EL constant in V
C1_EL = -0.1765;            % EL constant in V/C
C2_EL = 5.5015;             % EL constant in V
I_EL0 = 0.1341;             % EL constant in A
R_EL = -3.3189;             % EL constant in ohm.C
C_H2 = 2.39;                % conversion factor in Ah/L
n_I_EL = 0.7;               % EL efficiency
% in Kelouwani model hydrogen pressure is considered to be 1 atm, while in
% fact we use compressures, the energy conversion stayes true.
dH = 286;                   % Conversion constant (kJ/mol)
V_T = 22.4;                 % Conversion constant (L/mol)

% max. current limit of 24  PEM-EL cell
if I_EL/N_EL>(5e3/n_I_EL/Vbus/24)
    error('Over discharge of EL');
    I_EL = (5e3/n_I_EL/Vbus)*fix(N_EL/24);
end

% Calculation
Vel = V_EL0*ones(1,T) + C1_EL*T_EL + C2_EL*log(I_EL/I_EL0) + R_EL*I_EL./T_EL;
math_error_indx = Vel==-inf;
Vel(math_error_indx) = 0;

HGR = N_EL*(n_I_EL/C_H2)*(I_EL*Vbus/Vel*n_inv);
HGR(math_error_indx) = 0;
P_el = HGR*dH/V_T/3600;          % EL power in kW
end