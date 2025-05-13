clear
% clc
close all
%% SC2 article
% Depot_demand = [808.931426041718	808.931426041718	808.931426041718	808.931426041718	808.931426041718	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	808.931426041718	808.931426041718	808.931426041718];
% demandT1 = [0	0	0	0	0	0	57.2638888888889	187.833333333333	146.611111111111	212.722222222222	151.958333333333	235.472222222222	133	318.500000000000	144.083333333333	236.638888888889	227.305555555556	149.527777777778	311.888888888889	162.361111111111	0	0	0	0];
% demandT2 = [0	0	0	0	0	0	36	97.8500000000000	88.6000000000000	82.1000000000000	157.600000000000	153.150000000000	82.5000000000000	135.500000000000	129.050000000000	77.4500000000000	129.050000000000	141.150000000000	89.3000000000000	106.900000000000	0	0	0	0];
% curi_demand_line1 = [13838	13364	13119	13073	13273	13551	14196	15252	15956	16469	16619	16280	16576	16833	16837	16710	16490	16999	17041	17039	16921	16433	15568	14548]; 
%% SC1 article
Depot_demand = [1112.56080342936	1112.56080342936	1112.56080342936	1112.56080342936	1112.56080342936	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1112.56080342936	1112.56080342936	1112.56080342936];
demandT1 = [0	0	0	0	0	0	40.9544444444445	51.6666666666667	106.433333333333	89.0388888888889	77.2933333333333	107.880000000000	123.827777777778	89.6588888888889	59.0722222222222	82.0466666666667	71.0244444444444	52.1833333333333	0	32.4811111111111	0	0	0	0];
demandT2 = [0	0	0	0	0	0	83.9066666666667	41.7466666666667	91.6222222222222	89.4866666666667	90.6922222222222	75.2611111111111	99.8200000000000	63.6188888888889	93.9988888888889	60.4500000000000	99.1655555555556	24.8000000000000	0	28.3822222222222	0	0	0	0];
curi_demand_line1 = [13838	13364	13119	13073	13273	13551	14196	15252	15956	16469	16619	16280	16576	16833	16837	16710	16490	16999	17041	17039	16921	16433	15568	14548]; 


% Depot_demand = [808.931426041718	808.931426041718	808.931426041718	808.931426041718	808.931426041718	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	808.931426041718	808.931426041718	808.931426041718];
% demandT1 = [0	0	0	0	0	0	57.2638888888889	187.833333333333	146.611111111111	212.722222222222	151.958333333333	235.472222222222	133	318.500000000000	144.083333333333	236.638888888889	227.305555555556	149.527777777778	311.888888888889	162.361111111111	0	0	0	0];
% demandT2 = [0	0	0	0	0	0	36	97.8500000000000	88.6000000000000	82.1000000000000	157.600000000000	153.150000000000	82.5000000000000	135.500000000000	129.050000000000	77.4500000000000	129.050000000000	141.150000000000	89.3000000000000	106.900000000000	0	0	0	0];
% curi_demand_line1 = [13838	13364	13119	13073	13273	13551	14196	15252	15956	16469	16619	16280	16576	16833	16837	16710	16490	16999	17041	17039	16921	16433	15568	14548]; 

Pl = demandT1 + demandT2 + Depot_demand + curi_demand_line1;

LB = -100000*ones(1,24);  % Lower bound of variables
UB = 350000*ones(1,28);
LB(25:28) = 0;

x = optimvar('x', 1, 28, 'LowerBound', LB, 'UpperBound', UB);
[total_cost, pgrid, p_pv, PotherGen, p_wt dccc_bat] = hazine2(x,Depot_demand,demandT1,demandT2,curi_demand_line1); 
Mu=[0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 1 1.6 1.6 1.6 1 0.6 0.6]*0.2;


constraint1 = x(28)<=0.927*x(25);
constraint2 = x(1:24)<=0.75*x(25)*0.971;
constraint3 = x(1:24)>=-0.75*x(25)/0.927;
cons = [];
for i=1:24
    const1(i) = x(28) + sum(x(1:i));
    const2(i) = x(28) + sum(x(1:i));
end

% Set Problem to minimum
prob = optimproblem('ObjectiveSense','min');

prob.Objective = total_cost;    % Add Costfuncion to Problem
prob.Constraints.cons1 = constraint1;  % Add Constraint1 to Problem
prob.Constraints.cons2 = constraint2;  % Add Constraint2 to Problem
prob.Constraints.cons3 = constraint3;  % Add Constraint2 to Problem
prob.Constraints.cons5 = const1>=0.05*x(25);  % Add Constraint2 to Problem
prob.Constraints.cons6 = const2<=0.95*x(25);  % Add Constraint2 to Problem
prob.Constraints.cons7 = sum(x(1:24)) == x(28);  % Add Constraint2 to Problem
prob.Constraints.cons8 = PotherGen <= max(curi_demand_line1);
prob.Constraints.cons9 = p_pv - x(1:24) + p_wt <= 1.1*max(Pl);
x0.x = UB;
options = optimoptions(@fmincon,'MaxIterations',15000);
[sol,fval] = solve(prob,x0,'Options',options);
X = sol.x;
SOC(1) = X(28) / X(25);
for i=1:24
    SOC(i+1) = SOC(i) + X(i)/ X(25);
end
Pbat = X(1:24);
[total_cost pgrid p_pv PotherGen p_wt dccc_bat] = hazine2(X,Depot_demand,demandT1,demandT2,curi_demand_line1); 
PotherGen(PotherGen<0) = 0;

constraint4 = SOC<0.95;
D2G = pgrid.*(pgrid<0);
PV2G = D2G - 0.97*Pbat.*(Pbat<0);
G2D = pgrid.*(pgrid>0);
bar(G2D)
hold on
bar(Depot_demand + demandT1 + demandT2)
hold on
bar(PV2G,'g')
hold on
plot(SOC*X(25),'r','LineWidth',2)
hold on
plot(p_pv,'g-.','LineWidth',2)
hold on
bar(PotherGen)
% plot(PotherGen,'LineWidth',2)
plot(curi_demand_line1,'LineWidth',2)
plot(curi_demand_line1 + Depot_demand + demandT1 + demandT2,'LineWidth',2)
yyaxis right
plot(Mu,'LineWidth',2)
legend('G2D','Bus Demand','D2G','Battery Energy','PV','Other Generation','City Demand','Total Load')

disp(['Total Profit = ',num2str(total_cost)])
disp(['---Const 1---'])
disp(['p_pv - x(1:24) + p_wt  <= Pl: --- ','Yes'])
disp(['---Const 2---'])
disp(['p_pv - x(1:24) + p_wt <= 1.1*max(Pl): --- ','Yes'])
disp(['--------Other Output------'])
disp(['Capacity of PV: ',num2str(X(27)), ' with Max of 100000'])
disp(['Capacity of Bat: ',num2str(X(25)), ' with Max of 100000'])
disp(['Capacity of Wind: ',num2str(X(26)), ' with Max of 100000'])
disp(['Max of Other generation: ',num2str(max(PotherGen))])
disp(['Max of City Demand: ',num2str(max(curi_demand_line1))])