function PV = P_PV(G,T,n_inv)
% PV = P_PV(G,T,n_inv)
% G is insolation is kW/m2
% T is air temperaure in C
% P_PV will be one panel output in kW

% PV properties

G = 1000 * G;     % irradiation in W/m2

NOCT = 45;        % nominal operating cell temperature in C
Voc = 21;         % open circuit volatge in V
Isc = 6.5;        % short circuit current in A
K = 1.38e-23;     % Boltzman constant in J/K
q = 1.6e-19;      % electron charge in C
G0 = 1000;        % rated insolation in W/m2
T0 = 298;         % ambient temperature in K
Rs = 0.012;       % series resistant in ohm
B = 0.058;        % PV module technology specific-related coefficient 
a = 1.21;         % non-linearity effect coefficient
y = 1.15;         % non-linearity temperature-voltge related effect coefficient
nMPPT = 0.95;     % ideality factor  
noth = 0.95;      % other losses factor (such as dust on the panels)
nMPP = 1.17;      % coefficient


% calculation

T = T + (NOCT - 20)/800.*G;     % calculating module temperature in C
warning('caution: NOCT operator is applied');
T = 273*ones(1,size(T,2)) + T;  % module temperature in K    

for i=1:size(G,2)

    C(i) = Voc ./ (nMPP * K .* T(i) ./ q); % dummy variable

    % generated power by one panel
    P_panel(i) = ( C(i) - log( C(i) + 0.72) ) ./ ( 1 + C(i) )...
        .* ( 1 - Rs/Voc*Isc ) .* Isc .* ((G(i)./G0).^a)...
        .* ( Voc ./ (1 + B .* log( G0/G(i) ) ) )...
        .* ( ( T0./T(i) ).^y ) .* nMPPT .* noth;

end

PV.Ppv = n_inv.*P_panel./1000;
PV.Ipv = PV.Ppv./48;                    % DC bus voltage is 48V

end

