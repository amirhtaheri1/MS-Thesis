function check = Constraints(x,data,detail)
warning('BYPASS: Constraints.m for temporary testing');
check = true;
return
% check = Constraints(x,data,detail)
% modified for EV and FCV
% x = particle(i).position /in EPSO_adapted.m
% checkes system technical constraints during the simulation period T.
% if true (check=1) the solution in a feasible one.

%% Extracting data
% Each unit output and load from main.m (function: EPSO_adapted.m)
% Ppv = data.PV.Ppv;     % power output of PV in kW
% Pwt = data.WT.Pwt;     % power output of WT in kW
% load = data.load;   % load vector in kW

%% Initialization
% data extraction
N_PV = x(1);
N_WT = x(2);
E_BT_MAX = x(3);
C_SC = x(4);
E_H2_MAX = x(5);
N_EL = x(6);
N_FC = x(7);

%% Size constraints
N_PV_max = data.x_max(1);
N_WT_max = data.x_max(2);
E_BT_MAX_max = data.x_max(3);
C_SC_max  = data.x_max(4);
E_H2_MAX_max = data.x_max(5);
N_EL_max = data.x_max(6);
N_FC_max = data.x_max(7);

N_PV_min = data.x_min(1);
N_WT_min = data.x_min(2);
E_BT_MAX_min = data.x_min(3);
C_SC_min  = data.x_min(4);
E_H2_MAX_min = data.x_min(5);
N_EL_min = data.x_min(6);
N_FC_min = data.x_min(7);

LPSP_max = 0.10;                % max loss of energy supply (%)
LHSP_max = 0.30;                % max loss of hydrogen supply (%)

%% Detail of each system operating
% detail = details(x,data);       % calling details/detail of each system

BEV = data.BEV_flag;
FCEV = data.FCEV_flag;

SOC_BT = detail.SOC_BT;
SOC_SC = detail.SOC_SC;
LPSP = detail.LPSP;
LHSP = detail.LHSP;

SOC_BEV = detail.SOC_BEV;
LES = detail.LES;
LHS = detail.LHS;

SOC_BT_min = 0;                 % min SOC of the battery based on its capacity
SOC_BEV_min = 15;               % min SOC of EV
LES_max = 0;                    % max loss of BEV charge in days
LHS_max = 0;                    % max loss of FCEV charge in days

%% Constraint check
if BEV
    if ((N_PV>N_PV_max) || (N_WT>N_WT_max) || (E_BT_MAX>E_BT_MAX_max) || ...
        (C_SC>C_SC_max) || (E_H2_MAX>E_H2_MAX_max) || (N_EL>N_EL_max) || ...
        (N_FC>N_FC_max) || ...
        (N_PV<N_PV_min) || (N_WT<N_WT_min) || (E_BT_MAX<E_BT_MAX_min) || ...
        (C_SC<C_SC_min) || (E_H2_MAX<E_H2_MAX_min) || (N_EL<N_EL_min) || ...
        (N_FC<N_FC_min) || ...
        (min(SOC_BT)<SOC_BT_min) || (LPSP>LPSP_max) || (max(SOC_SC)>100) ||...
        (LHSP>LHSP_max) || (LES>LES_max) || (min(SOC_BEV)<SOC_BEV_min))

        check = false;
    else
        check = true;
    end

elseif FCEV
    if ((N_PV>N_PV_max) || (N_WT>N_WT_max) || (E_BT_MAX>E_BT_MAX_max) || ...
        (C_SC>C_SC_max) || (E_H2_MAX>E_H2_MAX_max) || (N_EL>N_EL_max) || ...
        (N_FC>N_FC_max) || ...
        (N_PV<N_PV_min) || (N_WT<N_WT_min) || (E_BT_MAX<E_BT_MAX_min) || ...
        (C_SC<C_SC_min) || (E_H2_MAX<E_H2_MAX_min) || (N_EL<N_EL_min) || ...
        (N_FC<N_FC_min) || ...
        (min(SOC_BT)<SOC_BT_min) || (LPSP>LPSP_max) || (max(SOC_SC)>100) ||...
        (LHSP>LHSP_max) || (LHS>LHS_max))

        check = false;
    else
        check = true;
    end
end
    
%% Logging results (with total cost(TC))
%{
if data.excel==1
    filename = 'details_FLC_log.xlsx';
    sheet = 1;
    xlRange = 'A1';
    
    % loading the excel file to retrieve previous results
    last_file = xlsread(filename);
    
    % changing the first value of the last log(written by CostFunction.m)
    new_file = [last_file(1:end-1,:);check,last_file(end,2:end)];
    xlswrite(filename,new_file,sheet,xlRange);
end
%}
end

