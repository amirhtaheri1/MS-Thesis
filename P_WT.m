function WT = P_WT(v,n_inv)
% WT = P_WT(v,n_inv)
% v is wind speed in m/s.
% P_WT will be wind turbine output power in kW.
% Pwt will be wind turbine output in kW.
% I wt will be wind turbine current in A.

% wind in hub height
%{
hub_height = 50;        % in m
anemometer_height = 50; % in m
a = 0.14;               % Hellman exponent / friction coefficient
v = v * (hub_height/anemometer_height)^a;
%}

% WT properties    
Uc = 2.5;               % cut-in speed in m/s
Ur = 9.5;               % rated speed in m/s
Uf = 25;                % furling speed in m/s
Pmax = 10;              % rated power in kW
Power = zeros(1,length(v));

% finding power for each wicd speed
for i=1:length(v)

    if (v(i)<Uc) || (v(i)>Uf)
        Power(i) = 0;

    elseif (v(i)>=Ur) && (v(i)<=Uf)
        Power(i) = Pmax;

    elseif (v(i)>=Uc) && (v(i)<Ur)
        Power(i) = Pmax * (v(i)-Uc)/(Ur-Uc);
    end

end

WT.Pwt = n_inv*Power;
WT.Iwt = Power./48;           % DC bus voltage is 48V

end

