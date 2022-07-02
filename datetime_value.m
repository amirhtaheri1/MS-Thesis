function output = datetime_value()
% output = datetime_value()
% current PC time
% output will be the datetime value in a series number.
% example: 20220528120539 for 2022, May, 28th, 12:05:39 p.m..

c = fix(clock);

% numerical value of time
output = c(6) + c(5)*1e2 + c(4)*1e4 + c(3)*1e6 + c(2)*1e8 + c(1)*1e10;
end