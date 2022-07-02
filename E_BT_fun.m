function BT = E_BT_fun(E_bt0,I_bt,E_BT_max,n_BT,n_inv)
% BT = E_BT_fun(E_bt0,I_bt,E_BT_max,n_BT,n_inv)
% Called from details_FLC.m (version 12).
% E_bt0 is battery energy in Ah.
% n_BT is battery efficiency.
% I_BT is battery current in A.
% E_BT will be battery energy in Ah.
% SOC_BT will be battery state of charge in %.
% I_BT will be battery charging current (+: charging, -: dischrging).
% I_SC will be the current more than I_max of battery (I_balance-I_BT).

% Initiallizarion
Ebt = 0;
Ibt = 0;
Isc = 0;
I_max = 0.2*E_BT_max;

% max. battery energy limit
% For safety check. This condition should be false.
if E_bt0>E_BT_max
%     error('battery capacity error');
    Ibt = 0;
    Ebt = E_bt0;
else
    I_abs = abs(I_bt);
    % battery charging (+) or discharging (-)
    dir = sign(I_bt);
    %{
        % n is charging of discharging efficiency
        if dir==1
            n=n_inv;
        elseif dir==-1
            n=n_BT*n_inv;
        end
    %}
    n=1;
    % 20_percent_max BT charge/discharge current limit
    if I_abs>I_max
        Ibt = n*dir*I_max;
        % SC charge(+) or deploy (-)
        Isc = n_inv*dir*(I_abs-I_max);
        Ebt = E_bt0 + Ibt;
    else
        Ibt = n*dir*I_abs;
        Ebt = E_bt0 + Ibt;
    end
    
    if Ebt>E_BT_max
        Ibt = Ibt - (Ebt-E_BT_max*dir);
        Ebt = E_BT_max*dir;
    end
    
    if Ebt<0
        Ibt = Ibt-Ebt;
        Ebt = 0;
        %             fprintf(['over discharge in i=',num2str(i),' prevented.\n']);
    end
    
end

% results
BT.E_BT = Ebt;                          % BT energy in Ah
BT.SOC_BT = 100 * BT.E_BT ./ E_BT_max;  % BT SOC
BT.I_BT = Ibt;
BT.I_SC = Isc;
end