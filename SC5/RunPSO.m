clc; clear; close all; warning off all;

%% PSO parameters
load integrated_data2.mat

 
cnt = 1;
for i=1:size(StNumMat,2)
    for j=1:size(StNumMat{i},2)
        NumTerminalSeen(cnt) = numTermialSeen{i}{j};
        cnt = cnt + 1;
    end
end
ChargeSchVar = sum(NumTerminalSeen);
ChargeSchVar_LB = zeros(1,ChargeSchVar);
ChargeSchVar_UB = ones(1,ChargeSchVar);
Costfunc = @(x) Cost(x);
Charger_types = 4;
nvars = ChargeSchVar + 5 + 1; 
LB = [100 0 0 1 1 ChargeSchVar_LB 1];
UB = [800 7 7 Charger_types Charger_types ChargeSchVar_UB 4];

% LimitSOC = 5 + 15.6;
% UpperSOC = 95;
% 
% LB = [300 LimitSOC LimitSOC 2 2 LimitSOC LimitSOC LimitSOC LimitSOC 10000 25000]; 
% UB = [900 UpperSOC UpperSOC 7 7 UpperSOC UpperSOC UpperSOC UpperSOC 25000 45000];

options = optimoptions('particleswarm','TolFun',1e-7,'PlotFcn','pswplotbestf');
rng default
tic;
[x,fval] = particleswarm(Costfunc,nvars,LB,UB,options);
toc;

y = round(x);
Cost_Optim = fval;

%%
