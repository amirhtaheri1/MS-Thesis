function [J,detail] = CostFunction(x,data)
% [J,detail] = CostFunction(x,data)
% cost = CostFunction_MPC(X,U,e,data,Ts)
% cost is in $/sample_week.
% new J is set by the formulatiob on (1401.03.09).
% from (1401.03.05) to (1401.03.09) (LPSP) upgraded to (LPSP*168/5).
% x(1:7) is equipment sizes (ver. 17)

% disp('CostFunction running...')

% CostFunction is running for BEV
BEV_flag = data.BEV_flag;
% CostFunction is running for FCEV
FCEV_flag = data.FCEV_flag;

%% Initialization
% simulation time steps
T = data.T;

% equipment sizes (in vector format)
N_PV = x(1);
N_WT = x(2);
E_BT_MAX = x(3);
C_SC = x(4);
E_H2_MAX = x(5);
N_EL = x(6);
N_FC = x(7);

Vbus = 48;
% years to simulate cost for
n = 1/52; % 1 week is n years
% WT rated power
P_WT_r = 10;               % WT rated power in kW
% max power of BT for an hour
P_BT_max = E_BT_MAX*Vbus*0.2/1000;
% SC max power
V_SC_max = 24.5;              % voltage of SC when charged in V
E_SC_max = 0.5 .* C_SC .* V_SC_max.^2/3.6E6;   % max. SC energy in kWh
% simulation step of an hour would result in equal power and energy for SC
% currently is dt.
dt = 6.225; % SC application time in seconds
P_SC_max = E_SC_max/dt;

%% life time for each system
n_WT = 25;          % wind turbine life time in yr
n_PV = 20;          % photovoltaic panel life time in yr
n_BT = 5;           % battery life time in yr
n_SC = 20;          % supercapacitor life time in yr
n_EL = 20;          % electrolyzer life time in yr
n_FC = 10;          % fuel cell life time in yr
n_HT = 15;          % hydrogen tank life time in yr
n_conv_inv = 20;    % converter life time in yr
n_comp = 20;        % compressor life time in yr

%% Capital costs
Ccap_WT = 2500;             % $/kW
Ccap_PV = 1000;             % $/kW
Ccap_BT = 200;              % $/kWh
Ccap_SC = 300;              % $/kWh
Ccap_EL = 2000;             % $/kW
Ccap_FC = 1500;             % $/kW
Ccap_HT = 30;               % $/kWh
Ccap_conv_inv = 200;        % $/kW
Ccap_comp = 46600;          % $

%% loan data
n_week = 1/52;      % simulation time in yr (1 week = 1/52 yr)
n = 20;             % loan payback time in yr
i = 0.05;           % loan interest

%% CRF and SPPF calculations
CRF = i * ( (1+i)^n ) / ( (1+i)^n -1); % capital recovery factor

k1 = n_BT : n;      % years needing to change batteries after
k2 = n_FC : n;      % years needing to change fuel cells after
k3 = n_HT : n;      % years needing to change hydrogen tanks after

SPPF1 = (1 + sum(1./(1+i).^k1));  % single payment present factor for BT
SPPF2 = (1 + sum(1./(1+i).^k2));  % single payment present factor for FC
SPPF3 = (1 + sum(1./(1+i).^k3));  % single payment present factor for HT

Cmod_BT = Ccap_BT * SPPF1;        % BT price with respect to SPPF
Cmod_FC = Ccap_FC * SPPF2;        % FC price with respect to SPPF
Cmod_HT = Ccap_HT * SPPF3;        % HT price with respect to SPPF

%% O&M costs
% fixed cost
Com_f_PV = 9.52;             % $/kW/yr
Com_f_WT = 10;               % $/kW/yr
Com_f_BT = 5;                % $/kW/yr
Com_f_SC = 5;                % $/kW/yr
% variable cost
Com_v_FC = 0.02;             % $/kWh
Com_v_EL = 0.0045;           % $/kWh
Com_v_BT = 0.05;             % $/kWh
Com_v_SC = 0.005;            % $/kWh
% Hydrogen equipment!, battery and converters have negligible O&M costs
% at least currently

%% Cost function
% calling details to get detail of each system
detail = details(x,data);

% capital and replacement cost
Ccap = CRF*(N_PV*Ccap_PV*max(data.PV.Ppv) + N_WT*P_WT_r*Ccap_WT + ...
            E_BT_MAX*Cmod_BT*Vbus/1000 + E_SC_max*Ccap_SC + ...
            E_H2_MAX*Cmod_HT + data.P_max_EL*Ccap_EL + ...
            data.P_max_FC*Cmod_FC + detail.P_Conv_inv*Ccap_conv_inv + ...
            Ccap_comp);
% retracting high-pressure compressor for BEV (making it half price)
if data.BEV_flag
   Ccap = Ccap - CRF*Ccap_comp/2; 
end
% max rated Ccap_n to nomilalize Ccap
Ccap_n = 23500;    % $

% O&M cost
Com_f = n_week * (N_PV*Com_f_PV*max(data.PV.Ppv) + N_WT*P_WT_r*Com_f_WT +...
             P_BT_max*Com_f_BT + P_SC_max*Com_f_SC);

% E_FC is negative
Com_v = -sum(detail.E_FC)*Com_v_FC + sum(detail.E_EL)*Com_v_EL +...
        sum(abs(diff(detail.E_BT)))*Vbus/1000*Com_v_BT + ...
        sum(abs(diff(detail.E_SC)))*Com_v_SC;

% max rated Com to nomilalize Com
Com_n = 22;                          % $/week

% sample weekly cost
% J = w*(Com_f + Com_v)/Com_n + (1-w)*detail.LPSP;
% OMcost = Com_f + Com_v;
   
if data.OS
    W0 = 0.3;
    W1 = 0.1;
    W2 = 0.3;
    W3 = 0.1;
    W4 = 0.05;
    W5 = 0.1;
    W6 = 0.05;
    W7 = 0;

    J = W0*(Ccap/Ccap_n)+ W1*(Com_f + Com_v)/Com_n + W2*(detail.LPSP)*168/5 +...
        W3*(1 - sum(detail.HL)/100/T) + W4*(sum(abs(diff(detail.HL))))/100/T +...
        W5*(1 - sum(detail.SOC_BT)/100/T) + W6*(sum(abs(diff(detail.SOC_BT))))/100/T + ...
        W7*(sum(detail.E_dump)/(sum(data.PV.Ppv+data.WT.Pwt-data.load/data.n_inv)));
else
    W1 = 0.1;
    W2 = 0.3;
    W3 = 0.2;
    W4 = 0.1;
    W5 = 0.2;
    W6 = 0.1;

    J = W1*(Com_f + Com_v)/Com_n + W2*(detail.LPSP)*168/5 +...
        W3*(1 - sum(detail.HL)/100/T) + W4*(sum(abs(diff(detail.HL))))/100/T +...
        W5*(1 - sum(detail.SOC_BT)/100/T) + W6*(sum(abs(diff(detail.SOC_BT))))/100/T;
end
 
format
C_cap = Ccap/Ccap_n 
C_OM = (Com_f + Com_v)/Com_n
LPSP = (detail.LPSP)*168/5
C_HTE = (1 - sum(detail.HL)/100/T)
C_HTD = (sum(abs(diff(detail.HL))))/100/T
C_BTE = (1 - sum(detail.SOC_BT)/100/T)
C_BTD = (sum(abs(diff(detail.SOC_BT))))/100/T

%%
% disp('CostFunction ended.')
end

