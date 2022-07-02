clc
close all
clear all

disp('MegaRun_over_Main running...');
% tMegaStart = tic;

%% common signals
% optimal sizing step flag
MegaRunSignal.OS =1;
MegaRunSignal.MaxRunNumber = 1;
% linear relation control flag
MegaRunSignal.LRC_flag = 0;
% fuzzy logic control flag
MegaRunSignal.FLC_flag = 1;
% using initial guess in the heuristic optimization algorithms flag
MegaRunSignal.initial_guess_flag = 1;
% using dynamic weights based on R. Lorestani and M. Ardehali (E-PSO)
% (https://doi.org/10.1016/j.renene.2017.12.037)
MegaRunSignal.dynamic_weights_flag = 0;
% last file MaxRunNumber recorded in the root folder
MegaRunSignal.MegaRunNumebr = 200;
% max. try in running each program for testing
MaxTry = 1;
% time for the cpu to cool down between simulations in seconds
CoolDownTime = 200;

%% FCEV-PSO
for k=1:MaxTry
MegaRunSignal.MegaRunNumebr = MegaRunSignal.MegaRunNumebr +1;
MegaRunSignal.EPSO_flag = 0;
MegaRunSignal.BEV_flag = 0;

Main(MegaRunSignal)
pause(CoolDownTime)
end

%% BEV-PSO
for i=1:MaxTry
MegaRunSignal.MegaRunNumebr = MegaRunSignal.MegaRunNumebr +1;
MegaRunSignal.EPSO_flag = 0;
MegaRunSignal.BEV_flag = 1;

Main(MegaRunSignal)
pause(CoolDownTime)
end

%% FCEV-EPSO
for z=1:MaxTry
MegaRunSignal.MegaRunNumebr = MegaRunSignal.MegaRunNumebr +1;
MegaRunSignal.EPSO_flag = 1;
MegaRunSignal.BEV_flag = 0;

Main(MegaRunSignal)
pause(CoolDownTime)
end

%% BEV-EPSO
for j=1:MaxTry
MegaRunSignal.MegaRunNumebr = MegaRunSignal.MegaRunNumebr +1;
MegaRunSignal.EPSO_flag = 1;
MegaRunSignal.BEV_flag = 1;

Main(MegaRunSignal)
pause(CoolDownTime)
end

%% 
% tMegaEnd = toc(tMegaStart);

%% play sound when done
sound_duration = 1; % in minutes
sample_rate = 200e3; % samples per second
sound_length = sound_duration*60*sample_rate; % 10 min for (1000 sample/second) play rate
sound([rand(sound_length,1),rand(sound_length,1)],sample_rate);
sound([rand(0.9*sound_length,1),rand(0.9*sound_length,1)],0.9*sample_rate);

% tMegaEnd = toc(tMegaStart);
disp('MegaRun_over_Main ended!!!');