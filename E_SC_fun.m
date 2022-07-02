function SC = E_SC_fun(Vinitial,C,I_s,n_inv)
% SC = E_SC_fun(Vinitial,C,I_s,n_inv)
% SC is ideal (n_SC=1).
% I_s is receiving current in supercapacitor (SC) in A.
% C is equivalent capacity of SC in F.
% V_sc is voltage of SC before charge in V.
% Initial voltage of the SC is expected to be 48 V.
% Vsc will be voltage of SC after an hour of charge in V.
% Esc will be energy of SC after an hour of charge in kWh.
% I will be the current used by SC A.
% I_unused will be the current not used (from battery originally).

%% Initialization
Vfinal = 0;
I = 0;
Vbus = 48;
V_max = 24.5;                         % voltage of SC when charged in V
% down-voltage to capacitor rated voltage causes up-current
I_s = Vbus/V_max*I_s;
E_max = 0.5 .* C .* V_max.^2/3.6E6;   % max. SC energy in kWh
% charge/discharge duration inwhich SC is applied
dt = 6.225;

%% Main loop of calculation
Vfinal = Vinitial - dt*I_s/C;
%  SC limits
if Vfinal>V_max
    Vfinal = V_max;
elseif Vfinal<0
%     warning('SC overdischarged');
    Vfinal = 0;
end
I = C*(Vfinal-Vinitial)/dt;

%% Results
SC.Vsc = Vfinal;                  % returns SC voltage in V
Esc = n_inv.*0.5.* C.*Vfinal.^2;  % SC energy after an hour of charge in J
SC.E_SC = Esc/3.6E6;              % SC energy after an hour of charge in kWh
SC.SOC_SC = 100*SC.E_SC/E_max;
SC.I_SC = I;
% SC.I_unused = abs(I_s) - abs(I);

end