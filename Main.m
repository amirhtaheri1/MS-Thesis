function Main(MegaRunSignal) 
% clc
% clear all
% close all

% main body of code
% finding minimum of cost function (Com+LPSP)
% version 15: parallel PSOs + parameters fixed (CostFunction.m)
% version 16: R. Zahedi paper digram + new SC function (6.225s charge/discharge)
% version 17: parallel EPSO + BEV/FCEV charging model
% PV, WT, EL, FC, HT sizing is set according to zahedi et al.
% version 17: MPC is addded. but it is supposed to work independantly from
% details.m to avoid overcalculating teh overlap of the two function.
% version 18: update PSO as parPSO based on parEPSO in ver. 17, plus Fixing
% MPCoptimized.
% version 18: CostFunction.m is updated based on penalty cost method.
% version 18: constraints.m is bypassed, not effecting PSO anymore.
% version 19: Optimal sizing try based on optimal FLC from ver. 18.
% version 19.1: FC converter size within P_conv_inv got fixed in details.m
% (line 350).
% version 19.1: MPC is deleted.
% version 19: uploaded to github.
disp('Main running...');

% optimal sizing flag for CostFunction.m
data.OS = MegaRunSignal.OS;

%% File handling flags
        % writing detail_run number in a excel file to be globally accessible
        % data.detail_run_filename = 'detail_run_number.xlsx';
        % xlswrite('detail_run_number.xlsx',0,1,'A1');
        
% xlswrite('details_FLC_log.xlsx',zeros(1,134),1,'A1');
% delete 'Results.xlsx'
run_number = 1;      % initial value to track runs of this program
stuck_number = 0;    % initial value to keep track of bad initail guesses
% max_run = 1;         % the number of times to run the program
data.excel = 1;      % flag for printing answers in excel
% DATa = load('handel.mat');    % loading ending alarm sound 

%% Which vehicle test is running
data.BEV_flag = MegaRunSignal.BEV_flag;
data.FCEV_flag = 1 - data.BEV_flag;

%% Which controller is running
data.LRC = MegaRunSignal.LRC_flag;
data.FLC = MegaRunSignal.FLC_flag;

%% Extracting data from the excel file
% reading data from excel file
Table = xlsread("test_data.xlsx");

load = Table(2,:);          % extracting second row of table(load data)
wind_speed = Table(3,:);    % extracting third row of table(wind data)
temperature = Table(4,:);   % extracting forth row of table(temp. data)
irradiation = Table(5,:);   % extracting fifth row of table(sun data)

load_factor = 1;
load = load_factor * load;

%% PV and WT output
% system efficiencies are set in (details_FLC.m)
% inv/conv efficiency based on details.m
data.n_inv = 0.95;                         % inverter/converter efficiency
n_inv = data.n_inv;
PV = P_PV(irradiation,temperature,n_inv);  % power output of PV in kW
WT = P_WT(wind_speed,n_inv);               % power output of WT in kW

% packing output and load data to use in CostFunction.m/Constraints.m
% R. Zahedi PV and WT diagram
% (https://doi.org/10.1016/j.energy.2020.117935)
data.PV.Ppv = xlsread('PV_test.xlsx');
data.WT.Pwt = xlsread('WT_test.xlsx');
warning('PV and WT are derived from the paper, not the initial data');
data.load = load;
data.temperature= temperature;
data.T = length(temperature);
data.T = 24;
warning('testing: shortened T for testing to 24.');

%% Problem definition
%{
decision variables are [Npv(ea) Nwt(ea) E_BT-MAX(Ah) C_SC(F) E_H2_MAX(kWh)
N_EL(ea) N_FL(ea)
%}
problem.load = load;                   % laod vector in kW
problem.x_nVar = 7;                    % search space dimension of x
problem.param_nVar = 78;               % search space dimension of FLC parameters
problem.rule_nVar = 27;                % search space dimension of FLC rules

% lower bound of decision variables
problem.x_min = [0.5 0.5 600 10 1*33.2 10 15];
% upper bound of decision variables
problem.x_max = [5 5 6000 800 5*33.2 65 70];

problem.parameters_min = [-1*ones(1,7) zeros(1,10) -1*ones(1,7) zeros(1,10) zeros(1,10) -1*ones(1,7) zeros(1,10) zeros(1,10) -1*ones(1,7)];
problem.parameters_max = [ones(1,78)];

problem.rules_min = [ones(1,27)];
problem.rules_max = [3*ones(1,27)];

% packing parameter limits to use in CostFunction.m/Constraints.m
data.x_min = problem.x_min;
data.x_max = problem.x_max;
% % Com ratio to J in CostFunction.m
data.w = 0;

%% data used in details.m
data.Vbus = 48; % v
% efficiencies (moved from details.m)
data.n_I_EL = 0.7;          % electrolyzer efficiency
data.n_I_FC = 0.45;         % fuel cell efficiency
data.n_comp = 0.73;         % compressor efficiency/Humbolt State University
data.n_inv = data.n_inv;    % inverter/converter efficiency
data.n_BT = 0.85;           % battery/electrical energy storage efficiency
data.s = 0.0002;            % hourly self-discharge rate of battery
data.n_EV = 0.85;           % tesla model 3 battery charger effciency (drivingelectic.com)
% hydrogen energy
data.hydrogen_energy = 33.2; % kWh per 
% FCV charging percent needed in simulation in every night
data.VHL = 0.4; % vehicle hydrogen load in kg for each day of the week
% EV energy load (need) for a day
data.VEL = 7.8; % vehicle energy load in kWh for each day of the week 
% max. power of EL per 24 cells in kW
data.P_max_EL = 5;%/n_I_EL;
% max. power of FC per 25 cells in kW
% FC power and current is negative since it shows output current (power)
data.P_max_FC = 7;%;/n_I_FC;

%% PSO Parameters
PSO_parameters.stuck_step = 2;       % max number continual in unchanged J
PSO_parameters.try_it = 2;           % max number of repeat itteration
PSO_parameters.max_iteration = 4;    % max number of iterations

% nPop and Kw are found by sweeping in ver6 of this program
nPop = 4;                            % initial value for nPop
Kw1 = 0.4;                           % initial value for Kw for x
Kw2 = 0.4;                           % initial value for Kw for parameters
Kw3 = 0.4;                           % initial value for Kw for rules

PSO_parameters.nPop = nPop;          % population size

PSO_parameters.w_min = 0.3;          % min inertia coeficient
PSO_parameters.w_max = 1;            % max inertia coeficient
PSO_parameters.Kw1 = Kw1;            % shape-shift coefficient for w for x
PSO_parameters.Kw2 = Kw2;            % shape-shift coefficient for w for parameters
PSO_parameters.Kw3 = Kw3;            % shape-shift coefficient for w for rules

PSO_parameters.c1_min = 0.2;         % min personal acceleration coefficient
PSO_parameters.c1_max = 1;           % max personal acceleration coefficient
PSO_parameters.Kc1 = 0.4;            % shape-shift coefficient for c1

PSO_parameters.c2_min = 0.2;         % min social acceleration coefficient
PSO_parameters.c2_max = 1;           % max social acceleration coefficient
PSO_parameters.Kc2 = 0.4;            % shape-shift coefficient for c2

PSO_parameters.CR = 0.8;             % crossover rate (probability)
PSO_parameters.F_min = 0;            % minimum random factor
PSO_parameters.F_max = 1;            % maximum random factor

% flag for showing best solution in each iteration
PSO_parameters.show_iteration = false;

% using the un-optimized data (parameters+rule-set) 
% as initial data (particle(1)) in both EPSO and PSO algorithms
PSO_parameters.use_initial_guess = MegaRunSignal.initial_guess_flag;

% changing weights each iteration
PSO_parameters.dynamic_weights = MegaRunSignal.dynamic_weights_flag;
if not(MegaRunSignal.dynamic_weights_flag)
%     % OFLC
%     PSO_parameters.Kw1 = 0.2;
%     PSO_parameters.Kc1 = 0.3;
%     PSO_parameters.Kc2 = 0.6;
    % OS
    PSO_parameters.Kw1 = 0.1;
    PSO_parameters.Kc1 = 0.4;
    PSO_parameters.Kc2 = 0.4;
end

%% Initial FLC parameters and rules
% modified initial data based on Zahedi et al. - optimized results
% modification: P_Hydrogen (P) from un-optimized from R, Zahedi

data.rule_set_initial = [3 3 3 3 3 2 1 1 1,...
                         1 2 2 2 3 3 3 3 3,...
                         1 2 2 1 2 2 1 1 1];
data.parameters_initial = [-0.45232 0.02970 -0.88080 -.59459 -0.30264 0.02525 0.36817,...
                           0.15932 0.22178 0.36778 0.16155 0.23664 0.32761 0.72759 0.29926 0.49877 0.59856,...
                           -0.54212 -0.11271 -0.04051 0.64179 0.65458 0.18855 0.53990,...
                           0 0.01116 0.02427 0.00843 0.02813 0.02846 0.08277, 0.05540, 0.08335, 0.09550,...
                           0.05459, 0.11511, 0.31224, 0.09608, 0.37701, 0.65326, 0.85336, 0.46682, 0.51584, 0.89899,...
                           -0.85363, -0.14565, -0.56789, 0.28279, 0.39252, -0.39779, 1,...
                           0.00659, 0.0300, 0.06514, 0.05326, 0.06311, 0.07622, 0.09170, 0.04511, 0.05164, 0.07592,...
                           0.01142, 0.14680, 0.29076, 0.11759, 0.43269, 0.51114, 0.88726, 0.00583, 0.21351, 0.38861,...
                           -0.82093, -0.27445, -0.79571, 0.21790, 0.72071, 0.23152, 0.73224];


% % importing optimal fuzzy logic control for optimizing size from a spread
% sheet
% if data.OS
%     if data.BEV_flag
%         OFLC = xlsread('BEV_OFLC.xlsx');
%     else
%         OFLC = xlsread('FCEV_OFLC.xlsx');
%     end
%     data.parameters_initial = OFLC(1:78);
%     data.rule_set_initial = OFLC(79:105);
% end

%% Main loop
t_step = tic;

%% Calling optimization algorithm
% decision variables are [Npv(ea) Nwt(ea) E_BT-MAX(Ah) C_SC(F)
% E_H2_MAX(kWh) N_EL(ea) N_FL(ea)

size_init = [1, 1, 1200, 100, 100, 24, 35];
if data.OS
    size_init = [1.2, 1.4, 1200, 100, 100, 24, 35];
    warning('size_init is rigged as initial size for OS.') 
end

out.stuck = 1;
out.best_solution.cost = inf;

% Main while-loop
while (out.stuck==1) && (run_number<=MegaRunSignal.MaxRunNumber) && (data.FLC)
    % optimizaing the controller
    if MegaRunSignal.EPSO_flag
        out = parEPSO_OS(problem, PSO_parameters, data, size_init);
%         out = parEPSO_OP(problem, PSO_parameters, data, size_init);
    else
        out = parPSO_OS(problem, PSO_parameters, data, size_init);
%         out = parPSO_OP(problem, PSO_parameters, data, size_init);
    end
    
    % logging the timing
    test_time = toc(t_step)
    run_number = run_number +1;
    run_time = xlsread('run_time.xlsx');
    xlswrite('run_time.xlsx',[run_time;run_number,test_time,datetime_value],1,'A1');
end

% warning('Decision variables are [Npv(ea) Nwt(ea) E_BT-MAX(Ah) C_SC(F) E_H2_MAX(kWh) N_EL(ea) N_FL(ea)]')
warning('Results are saved in details_FLC_log.xlsx.')
% disp('Details_FLC.xlsx log is [stuck flag, w, bunch of zeros as spacer, optimal parameters, optimal rules]');
if out.stuck
    disp('Unvalid results!')
end
% end % starts at line 128 (for w=0.2:0.2:1)
% amid tun time
tEnd = toc(t_step)

% recover the break-it overflow problem in PSO/EPSO problem
if out.break_it==0
    out.break_it = PSO_parameters.max_iteration;
end

% fixing best_cost vector for plotting
out.best_costs = [out.best_costs(1:out.break_it); out.best_solution.cost*ones(PSO_parameters.max_iteration-out.break_it,1)];

% plotting best costs
figure; plot(out.best_costs,'k','LineWidth',1.05);
ylabel('Cost function J ($)'); xlabel('Iteration'); ytickformat('%0.3f');
ax = gca; ax.YAxis.Exponent = 0;

%% play sound when done
sound_duration = 0.5; % in minutes
sample_rate = 1e3; % samples per second
sound_length = sound_duration*60*sample_rate; % 10 min for (1000 sample/second) play rate
sound([rand(sound_length,1),rand(sound_length,1)],sample_rate);
sound([rand(1.2*sound_length,1),rand(1.2*sound_length,1)],1.2*sample_rate);
sound([rand(2*sound_length,1),rand(2*sound_length,1)],2*sample_rate);

%% saving .mat files for later-on reviews
% n1 = 'MegaRun_FLC_';
n1 = 'MegaRun_OS_';
if data.BEV_flag; n2='BEV_'; else n2='FCEV_'; end
if MegaRunSignal.EPSO_flag; n3='EPSO_'; else n3='PSO_'; end
if MegaRunSignal.initial_guess_flag; n4='initial-data_'; else n4='no-initial-data_'; end
if MegaRunSignal.dynamic_weights_flag; n5='dynamic-weights_'; else n5='fixed-weights_'; end
n6 = ['MegaRun_',num2str(MegaRunSignal.MegaRunNumebr),'_'];
n7 = ['Main-run-number_',num2str(run_number)];
MatFileName = strcat(n1,n2,n3,n4,n5,n6,n7);

save(MatFileName);

%% return
disp('Main ended.')
end