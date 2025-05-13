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

options = optimoptions('particleswarm','TolFun',1e-7,'UseParallel',true,'PlotFcn','pswplotbestf');
rng default
tic;
[x,fval] = particleswarm(Costfunc,nvars,LB,UB,options);
toc;

y = round(x);
Cost_Optim = fval;
