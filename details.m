function detail = details(x,data)
% detail = details(x,data)
% version 16: R. Zahedi paper digram + new SC function (6.225s charge/discharge)
% version 14, P_max_FC=7 kw istead of 5kW in my MS project.
% version 13, code (+equipment functions called) is renewed for series 
% caculation.
% version 12, code is renewed based on FLC of R.Zahedi and M.Ardehali
% model.
% version 19.1: max(E_FC) got corrected to abs(min(E_FC)) in P_conv_inv
% (line 350);
% EV: Tesla Model 3 (2021), 80 kwh battery, 26 kWh/100mile efficiency.
% FCEV: Toyota Mirai (2021), 5.6 kg H2_tank, 45.5 kWh/100mile efficiency.
% Based on x (particle(i).position) from EPSO_adapted and data from main.m.
% Simulation for 168 hour of a sample week.

%% Global counting of runs in detail_FLC
% detail_run = xlsread(data.detail_run_filename);
% detail_run = detail_run +1;
% xlswrite(data.detail_run_filename,detail_run,1,'A1');
% fprintf('run #%d of details_FLC.m\n',detail_run);

%% Size check
if length(x)~=(7+27+78)
    error('MAIN VECTOR SIZE IS NOT (7+27+78)!');
end

% disp('details_FLC running...')

%% Extracting data
% each unit output and load from main.m (function EPSO)
parameters = x(8:85);
rule_set = x(86:112);
x = x(1:7);

PV = data.PV;                       % power output of PV in kW
WT = data.WT;                       % power output of WT in kW
load = data.load;                   % load vector in kW
temperature = data.temperature;     % temperature in C
    
%% Initialization
% Simulation duration in hours
T = data.T;
% Sizes (in vector format)
N_PV = x(1);
N_WT = x(2);
E_BT_MAX = x(3);
C_SC = x(4);
E_H2_MAX = x(5);
N_EL = x(6);
N_FC = x(7); 

% updating max. power of EL and FC based on the number of their panels.
% originally were EL: 24 panels=5 kW and FC: 35 panels=7kW
data.P_max_EL = N_EL/24*5;      % in kW
data.P_max_FC = N_FC/35*7;      % in kW

% controller is running for EV
BEV_flag = data.BEV_flag;
% controller is running for FCV
FCEV_flag = data.FCEV_flag;
% initial value / loss of power supply
% LPS = 0;
% initial value / loss of hydrogen supply
LHS = 0;
% DC bus voltage in V
Vbus = data.Vbus;
% initial value / loss of EV battery supply (EV is not charged)
LES = 0;
% BT is initially full.
% SOC initial stage(100% charged)
SOC_BT = [];
SOC_BT(1) = 25; % initial charge of battery
E_BT = [];
E_BT(1) = SOC_BT(1)/100*E_BT_MAX;
% initial charge of EV battery
SOC_EV = [];
SOC_EV(1) = SOC_BT(1);
% initial SC abd BT state
I_BT = [];
I_BT_only = [];
I_SC = [];
Vinitial_SC = [];
% SC is initially empty (V_max=24.5V when its fully charged)
Vinitial_SC(1) = 0;
SOC_SC = [];
SOC_SC(1) = 0;

% SC energy
E_SC = [];
% FLC initial output
P_Hydrogen_FLC1 = [];
P_SC_FLC2 = [];
P_SC_FLC3 = [];
I_EL = zeros(1,T);
I_FC = zeros(1,T);
% Battery charge power initialization
E_EES = [];
% EV battery charge/discharge power
E_EV = E_EES;
% hydrogen energy
hydrogen_energy = data.hydrogen_energy; % kWh per kg
% HL initial stage(0% charged)
HL = [];
HL(1) = 50; % initial charge of hydrogen tank
% H2 tank energy
E_H2 = [];
E_H2(1) = HL(1)/100*E_H2_MAX;
% EV discharge power
EV = [];
% Electrolyzer power initialization
E_EL = zeros(1,T);
% fuel cell power initialization
E_FC = zeros(1,T);
% Power dumping initialization
E_dump = zeros(1,T);
% Power shortage initialization
E_shortage = zeros(1,T);
% FCV charging percent needed in simulation in every night
VHL = data.VHL; % vehicle hydrogen load in kg for each day of the week
% EV energy load (need) for a day
VEL = data.VEL; % vehicle energy load in kWh for each day of the week 
% EV battery capacity in kWh
C_EV = 80;
comp = 0;  % initial value for kg of hydrogen compressed
% FCEV charging power
FCEV = zeros(1,length(load));
% BEV charging power
BEV = zeros(1,length(load));
% FLC initial output valus
P_Hydrogen_FLC1 = [];
P_SC_FLC2 = [];
P_SC_FLC3 = [];

%% Efficiencies
n_I_EL = data.n_I_EL; % electrolyzer efficiency
n_I_FC = data.n_I_FC; % fuel cell efficiency
n_comp = data.n_comp; % compressor efficiency/Humbolt State University
n_inv = data.n_inv;   % inverter/converter efficiency
n_BT = data.n_BT;     % battery/electrical energy storage efficiency
s = data.s;           % hourly self-discharge rate of battery
n_EV = data.n_EV;     % tesla model 3 battery charger effciency (drivingelectic.com)

%% Max energy and power of EL, FC, and SC
% max. SC energy in kWh
E_SC_max = 0.5.*C_SC.*24.5.^2/3.6E6;  
% max. power of EL per 24 cells in kW
P_max_EL = data.P_max_EL;%/n_I_EL;
% max. power of FC per 25 cells in kW
% FC power and current is negative since it shows output current (power)
P_max_FC = data.P_max_FC;%;/n_I_FC;

%% Energy Balance
% Energy balance estimate of system in kW for the whole week
E_balance = N_PV*PV.Ppv + N_WT*WT.Pwt - load/n_inv;
Pn = max(E_balance);
dp = zeros(1,T);
I_BT_max = max(max(1000*E_balance/Vbus),-min(1000*E_balance/Vbus));
ctrl_input = zeros(T,4);

%% Main loop of calculation
% simulation for 168 hour of a sample week
for time=2:T
% % testing
% if (time==4)
%     keyboard
% end

% warning('currently using if-condition timing for charging, not P_EV vector')
% Simulating FCEV charging at 8 A.M. every day of the sample week
if (7<mod(time,24)) && (mod(time,24)<9) && FCEV_flag && (time>23)
    % if hydrogen supply suffices
    if (HL(time-1)/100*E_H2_MAX)>(hydrogen_energy*VHL)
        HL(time-1)=HL(time-1) - 100*(hydrogen_energy*VHL/E_H2_MAX);
        E_H2(time-1)=E_H2(time-1) - VHL*hydrogen_energy/n_comp;
        FCEV(time)=VHL*hydrogen_energy/n_comp;
    else
        LHS = LHS +1;
    end
end

% Simulating BEV charging at 12 A.M. until 6 P.M. every day of the sample week
if (0<mod(time,24)) && (mod(time,24)<7) && BEV_flag && (time>23)
    % if battery supply suffices. VEL is the energy for the 6 hours of
    % chrging which eqauls to VEL/6 gor each hour.
    if (SOC_BT(time-1)/100*E_BT_MAX*Vbus/1000)>(VEL/6)
        SOC_BT(time-1)=SOC_BT(time-1) - 100*VEL/6/E_BT_MAX/Vbus*1000;
        BEV(time)=VEL/6/n_inv;
        
    % special case: if SOC_BT is insufficient, but HT is suffcient
    % warning('special case needs to be tested!');
    elseif (HL(time-1)/100*E_H2_MAX)>(VEL/6)
        HL(time-1)=HL(time-1) - 100*(VEL/6/E_H2_MAX)/n_I_FC/n_inv;
        E_H2(time-1)=E_H2(time-1) - VEL/6/n_I_FC/n_inv;
        E_FC(time-1) = - VEL/6/n_I_FC/n_inv;
        BEV(time)=VEL/6/n_inv;
    else
        LES = LES +1;
    end
end
    
E_balance(time) = N_PV*PV.Ppv(time) + N_WT*WT.Pwt(time) - load(time)/n_inv;
dp(time) = E_balance(time)/Pn;  
% Vbus is consideed to be 48V.
I_BT_only(time) = 1000*E_balance(time)/Vbus;
% BT is initially full.
BT = E_BT_fun(E_BT(time-1),I_BT_only(time),E_BT_MAX,n_BT,n_inv);
E_BT(time) = BT.E_BT;
I_BT(time) = BT.I_BT;
% I_SC(time) = BT.I_SC/n_inv;
SOC_BT(time) = BT.SOC_BT;

% Normalizing I_BT
% I_BT_uni = I_BT(time)/I_BT_max;
I_BT_uni = I_BT_only(time)/I_BT_max;

% FLC and LRC controll input
ctrl_input(time,:) = [dp(time);SOC_BT(time-1)/100;I_BT_uni;SOC_SC(time-1)/100];

if data.FLC
    % fuzzy logic control (FLC)
    [P_Hydrogen_FLC1(time),P_SC_FLC2(time),P_SC_FLC3(time)] = FLC_optimized(ctrl_input(time,:),parameters,rule_set);
    % the sign direction of P_Hydrogen_FLC1 is applied like R.Zahedi.
end

if data.LRC
    % linear relation control (LRC)
    [P_Hydrogen_FLC1(time),P_SC_FLC2(time),P_SC_FLC3(time)] = LRC_unoptimized(ctrl_input(time,:));
end

% data preparation (I_EL, I_FC, and I_SC) in Amperes
% the sign direction of P_Hydrogen_FLC1 is applied like R.Zahedi.
% if P_Hydrogen_FLC1(time)>=0 produce electricity, else produce hydrogen.
if P_Hydrogen_FLC1(time)>=0
    I_FC(time) = -P_Hydrogen_FLC1(time)*P_max_EL/48*1000;
else
    I_EL(time) = -P_Hydrogen_FLC1(time)*P_max_FC/48*1000;
end
dt = 6.225/3600; % SC charge/discharge time in hour
I_SC_FLC2 = P_SC_FLC2(time)*E_SC_max*1000/48/dt;
I_SC_FLC3 = P_SC_FLC3(time)*E_SC_max*1000/48/dt;
 
if I_BT(time)>=0
    % FLC3 is deployed while BT is charging.
    I_SC(time) = I_SC_FLC3;
elseif I_BT(time)<0
    % FLC2 is deployed while BT is discharging.
    I_SC(time) = I_SC_FLC2;
end

% initiating SC for the next itteration
SC = E_SC_fun(Vinitial_SC(time-1),C_SC,I_SC(time),n_inv);
Vinitial_SC(time) = SC.Vsc;
SOC_SC(time) = SC.SOC_SC;
I_SC(time) = SC.I_SC;
E_SC(time) = SC.E_SC;

% H2 tank
H2 = E_H2_fun(temperature(time),I_EL(time),N_EL,n_I_EL,I_FC(time),N_FC,n_I_FC,E_H2_MAX,HL(time-1));
HL(time) = H2.HL;
E_H2(time) = H2.E_H2;
% energy used by EL (not produced!)
E_EL(time) = H2.E_EL;
% energy produced by FC
E_FC(time) = H2.E_FC;

% Energy balance after using BT, SC, EL, FC, H2
% BT efficiency is 1
% SC is divided by n_inv to encounter loss
% EL and FC energy contains loss
E_used(time) = Vbus/1000*(E_BT(time)-E_BT(time-1)) + (E_SC(time)-E_SC(time-1))/n_inv...
                + E_EL(time) + E_FC(time);
% E_FC and E_EL are produced and consumed energy, but E_BT, E_SC, and E_H2
% are energy available in these stor    age units
E_balance_after(time) = E_balance(time) - E_used(time);

% using EXTRA ENERGY in E_balance_after to charge battery
% moreover, use battery to fullfill lack of energy
I_BT_recheck(time) = 1000*E_balance_after(time)/Vbus;
BT2 = E_BT_fun(E_BT(time),I_BT_recheck(time),E_BT_MAX,n_BT,n_inv);
E_BT(time) = BT2.E_BT;
I_BT(time) = I_BT(time) + BT2.I_BT;
SOC_BT(time) = BT2.SOC_BT;

% updating E_balance_after
E_used(time) = Vbus/1000*(E_BT(time)-E_BT(time-1)) + (E_SC(time)-E_SC(time-1))/n_inv...
                + E_EL(time) + E_FC(time);
E_balance_after(time) = E_balance(time) - E_used(time);

% Matlab definite value error (Matlab zero is anything below 1e-13)
% consider values below 1e-9 as zero
if abs(E_balance_after(time))<1e-9
    E_balance_after(time) = 0;
elseif E_balance_after(time)>0
    %     disp('ENERGY DUMP NEEDED!');
    E_dump(time) = E_balance_after(time);
else
%     disp('details.m: ENERGY BALANCE VIOLATION! - LPS');
end

end

% Calculating reliability indices
LPS = E_balance_after<0;
LPSP = (LES+LHS+sum(LPS))/T;
LHSP = LHS/T/24;

% Power dump and shartage
% This line fo rE_dump gives the same results as the in-loop at the end of the Main loop
E_dump(~LPS) = E_balance_after(~LPS); 
E_shortage(LPS) = E_balance_after(LPS);
P_Conv_inv = sum([max(diff(E_SC/n_inv)),max(PV.Ppv/n_inv),...
                  max(WT.Pwt/n_inv),abs(min(E_FC)), max(E_EL/n_inv)]);
%                   max(WT.Pwt/n_inv),max(E_FC), max(E_EL/n_inv)]);
warning('in ver. 20 max(E_FC) is corrected into abs(min(E_FC)) compared to ver. 19');

% BT check
% negative_index = E_balance_after<0;
% if E_BT(negative_index)>0
%     ('did not fully use BT');
% end

%% Temperary results (no EV and FCEV)
% warning('parameters related to BEV & FCEV are constant and rigged in ver13');

% for Plotting_Results_from_log.m
detail.dp = dp;
detail.HL = HL;
detail.I_BT = I_BT;
detail.I_SC = I_SC;
detail.P_Hydrogen_FLC1 = P_Hydrogen_FLC1;
detail.P_SC_FLC2 = P_SC_FLC2;
detail.P_SC_FLC3 = P_SC_FLC3;

% for logging results in CostFunction.m
detail.P_max_EL = P_max_EL;
detail.P_max_FC = P_max_FC;
detail.E_balance = E_balance;
detail.E_balance_after = E_balance_after;
detail.I_BT_max = I_BT_max;

% for CostFunction.m
% Energy indices added in ver. 14 for variable O&M costs
detail.E_SC = E_SC;
detail.E_FC = E_FC;
detail.E_EL = E_EL;
detail.E_BT = E_BT;

% Power indices
detail.P_EL = max(diff(E_EL));
detail.P_FC = max(diff(E_FC));
detail.P_Conv_inv = P_Conv_inv;

% for Constraints.m
detail.SOC_BT = SOC_BT;
detail.SOC_SC = SOC_SC;
detail.LPSP = LPSP;
detail.LHSP = LHSP;
detail.SOC_BEV = 20*ones(1,T);
detail.LES = LES;
detail.LHS = LHS;

% for main.m
detail.E_dump = E_dump;
detail.E_shortage = E_shortage;
% detail

% EV charging
detail.BEV = BEV;
detail.FCEV = FCEV;

%%
% disp('details ended.')
end