function x = FixSize(x)
% x = FixSize(x)
% Fixing the equipment size of the HGPS used in parEPSO_OS.m and parPSO_OS.m.
% x is the size_init vector.
% The equipments are [PV, WT, BT, SC, HT, EL, FC].
% The equipment sizing will be the multiples of [0.01, 0.01, 10, 5,
% 6.66(0.2kg), 1, 1].
x(1) = 0.01*ceil(x(1)/0.01);
x(2) = 0.01*ceil(x(2)/0.01);
x(3) = 10*ceil(x(3)/10);
x(4) = 5*ceil(x(4)/5);
x(5) = 6.66*ceil(x(5)/6.66);
x(6) = 1*ceil(x(6)/1);
x(7) = 1*ceil(x(7)/1);
end
