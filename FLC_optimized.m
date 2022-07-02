function [output1,output2,output3] = FLC_optimized(input,parameters,rule_set)
% [output1,output2,output3] = FLC_optimized(input,parameters,rule_set)
% FLC for battery, Electrolyer, and Supercapacitor.
% Vector input is the FLC input values.
% Vector parameters is the parameters for all three FLCs.
% Vector rule_set is the rules for all three FLCs.
% FLC_unoptimized([dp(-1,1); SOC_BT(0,1); I_BT(-1,1); SOC_SC(0,1)])
% Output Will be P_Hydrogen_FLC1(-1,1),P_SC_FLC2(-1,1),and P_SC_FLC3(-1,1).
% Note: in FLC3 the input1mf is trapmf, but in FLC2 is trimf.

%% this part is intended for testing the FLC boundary calculation
if (logical(sum(input>=ones(4,1),'all'))||(input(2)<=0)||(input(4)<=0)||(input(1)<=-1)||(input(3)<=-1))
    % shifting boundaries so it contains (-1, 0, and 1)
    boundaries = input>=1;
    input(boundaries) = 0.99;
%     boundaries2 = input(1)<=-1;
%     input(boundaries2) = -0.99;
%     boundaries3 = input(2:4)<=0;
%     input(logical([0 boundaries3])) = 0.01;
    if input(1)<=-1; input(1) = -0.99; end
    if input(3)<=-1; input(3) = -0.99; end
    if input(2)<=0; input(2) = 0.01; end
    if input(4)<=0; input(4) = 0.01; end
end

%% Initiallization
% warning('off','all');
input1 = input(1:2)';
input2 = input(3:4)';

% FLC1
dp_p = parameters(1:7);
SOC_BT_p = parameters(8:17);
P_hydrogen_p = parameters(18:24);
rules1 = rule_set(1:9);

% FLC2
I_BT_p = parameters(25:34);
SOC_SC_p = parameters(35:44);
P_SC_p = parameters(45:51);
rules2 = rule_set(10:18);

% FLC3
I_BT_p2 = parameters(52:61);
SOC_SC_p2 = parameters(62:71);
P_SC_p2 = parameters(72:78);
rules3 = rule_set(19:27);

%% FLC1 for P_Hydrogen
fis = mamfis('Name',"FLC1");
fis = addInput(fis,[-1 1],'Name',"dp");

fis = addMF(fis,"dp","trimf",[-1 dp_p(1) dp_p(2)],"Name","N");
fis = addMF(fis,"dp","trimf",[dp_p(3) dp_p(4) dp_p(5)],"Name","Z");
fis = addMF(fis,"dp","trimf",[dp_p(6) dp_p(7) 1],"Name","P");

fis = addInput(fis,[0 1],'Name',"SOC_BT");
fis = addMF(fis,"SOC_BT","trapmf",[0 SOC_BT_p(1) SOC_BT_p(2) SOC_BT_p(3)],"Name","L");
fis = addMF(fis,"SOC_BT","trapmf",[SOC_BT_p(4) SOC_BT_p(5) SOC_BT_p(6) SOC_BT_p(7)],"Name","M");
fis = addMF(fis,"SOC_BT","trapmf",[SOC_BT_p(8) SOC_BT_p(9) SOC_BT_p(10) 1],"Name","H");

fis = addOutput(fis,[-1 1],'Name',"P_Hydrogen");
fis = addMF(fis,"P_Hydrogen","trimf",[-1 P_hydrogen_p(1) P_hydrogen_p(2)],"Name","N");
fis = addMF(fis,"P_Hydrogen","trimf",[P_hydrogen_p(3) P_hydrogen_p(4) P_hydrogen_p(5)],"Name","Z");
fis = addMF(fis,"P_Hydrogen","trimf",[P_hydrogen_p(6) P_hydrogen_p(7) 1],"Name","P");

% rulelist template: [input1mf input2mf outputmf ruleweight operator]
% 1 for the operator AND, 2 for the operator OR
rulelist = [1 1 rules1(1) 1 2
            1 2 rules1(2) 1 2
            1 3 rules1(3) 1 2
            2 1 rules1(4) 1 2
            2 2 rules1(5) 1 2
            2 3 rules1(6) 1 2
            3 1 rules1(7) 1 2
            3 2 rules1(8) 1 2
            3 3 rules1(9) 1 2];
fis = addRule(fis,rulelist);
output1 = evalfis(fis,input1);

if input(3)<0
% making I_BT a positive value.
input(3) = -input(3);
input2 = input(3:4)';

%% FLC2 for P_SC in BT discharging mode
fis2 = mamfis('Name',"FLC2");
fis2 = addInput(fis2,[0 1],'Name',"I_BT");

fis2 = addMF(fis2,"I_BT","trapmf",[0 I_BT_p(1) I_BT_p(2) I_BT_p(3)],"Name","N");
fis2 = addMF(fis2,"I_BT","trapmf",[I_BT_p(4) I_BT_p(5) I_BT_p(6) I_BT_p(7)],"Name","Z");
fis2 = addMF(fis2,"I_BT","trapmf",[I_BT_p(8) I_BT_p(9) I_BT_p(10) 1],"Name","P");

fis2 = addInput(fis2,[0 1],'Name',"SOC_SC");
fis2 = addMF(fis2,"SOC_SC","trapmf",[0 SOC_SC_p(1) SOC_SC_p(2) SOC_SC_p(3)],"Name","L");
fis2 = addMF(fis2,"SOC_SC","trapmf",[SOC_SC_p(4) SOC_SC_p(5) SOC_SC_p(6) SOC_SC_p(7)],"Name","M");
fis2 = addMF(fis2,"SOC_SC","trapmf",[SOC_SC_p(8) SOC_SC_p(9) SOC_SC_p(10) 1],"Name","H");

fis2 = addOutput(fis2,[-1 1],'Name',"P_SC");
fis2 = addMF(fis2,"P_SC","trimf",[-1 P_SC_p(1) P_SC_p(2)],"Name","N");
fis2 = addMF(fis2,"P_SC","trimf",[P_SC_p(3) P_SC_p(4) P_SC_p(5)],"Name","Z");
fis2 = addMF(fis2,"P_SC","trimf",[P_SC_p(6) P_SC_p(7) 1],"Name","P");

% rulelist template: [input1mf input2mf outputmf ruleweight operator]
% 1 for the operator AND, 2 for the operator OR
rulelist2 = [1 1 rules2(1) 1 2
             1 2 rules2(2) 1 2
             1 3 rules2(3) 1 2
             2 1 rules2(4) 1 2
             2 2 rules2(5) 1 2
             2 3 rules2(6) 1 2
             3 1 rules2(7) 1 2
             3 2 rules2(8) 1 2
             3 3 rules2(9) 1 2];
fis2 = addRule(fis2,rulelist2);
output2 = evalfis(fis2,input2);
output3 =0;

elseif input(3)>=0
%% FLC3 for P_SC for BT charging
fis3 = mamfis('Name',"FLC3");
fis3 = addInput(fis3,[0 1],'Name',"I_BT");

fis3 = addMF(fis3,"I_BT","trapmf",[0 I_BT_p2(1) I_BT_p2(2) I_BT_p2(3)],"Name","N");
fis3 = addMF(fis3,"I_BT","trapmf",[I_BT_p2(4) I_BT_p2(5) I_BT_p2(6) I_BT_p2(7)],"Name","Z");
fis3 = addMF(fis3,"I_BT","trapmf",[I_BT_p2(8) I_BT_p2(9) I_BT_p2(10) 1],"Name","P");

fis3 = addInput(fis3,[0 1],'Name',"SOC_SC");
fis3 = addMF(fis3,"SOC_SC","trapmf",[0 SOC_SC_p2(1) SOC_SC_p2(2) SOC_SC_p2(3)],"Name","L");
fis3 = addMF(fis3,"SOC_SC","trapmf",[SOC_SC_p2(4) SOC_SC_p2(5) SOC_SC_p2(6) SOC_SC_p2(7)],"Name","M");
fis3 = addMF(fis3,"SOC_SC","trapmf",[SOC_SC_p2(8) SOC_SC_p2(9) SOC_SC_p2(10) 1],"Name","H");

fis3 = addOutput(fis3,[-1 1],'Name',"P_SC");
fis3 = addMF(fis3,"P_SC","trimf",[-1 P_SC_p2(1) P_SC_p2(2)],"Name","N");
fis3 = addMF(fis3,"P_SC","trimf",[P_SC_p2(3) P_SC_p2(4) P_SC_p2(5)],"Name","Z");
fis3 = addMF(fis3,"P_SC","trimf",[P_SC_p2(6) P_SC_p2(7) 1],"Name","P");

% rulelist template: [input1mf input2mf outputmf ruleweight operator]
% 1 for the operator AND, 2 for the operator OR
rulelist3 = [1 1 rules3(1) 1 2
             1 2 rules3(2) 1 2
             1 3 rules3(3) 1 2
             2 1 rules3(4) 1 2
             2 2 rules3(5) 1 2
             2 3 rules3(6) 1 2
             3 1 rules3(7) 1 2
             3 2 rules3(8) 1 2
             3 3 rules3(9) 1 2];
fis3 = addRule(fis3,rulelist3);
output3 = evalfis(fis3,input2);
output2 =0;
end
%% Return
% warning('on','all');

end