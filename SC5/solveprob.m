clear
% clc
close all
Depot_demand = [268.105182200272	268.105182200272	268.105182200272	268.105182200272	268.105182200272	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	268.105182200272	268.105182200272	268.105182200272];
demandT1 = [0	0	0	0	0	0	0	90	269.250000000000	280	303	356.375000000000	397.500000000000	420.250000000000	199.250000000000	184.375000000000	286.625000000000	181.875000000000	449.375000000000	133.750000000000	0	0	0	0];
demandT2 = [0	0	0	0	0	0	0	90	386.375000000000	166.750000000000	339.625000000000	353.375000000000	320.875000000000	420.250000000000	326	378	434.750000000000	189.125000000000	449.250000000000	46.8750000000000	0	0	0	0];
curi_demand_line1 = [13838	13364	13119	13073	13273	13551	14196	15252	15956	16469	16619	16280	16576	16833	16837	16710	16490	16999	17041	17039	16921	16433	15568	14548]; 

Pl = demandT1 + demandT2 + Depot_demand + curi_demand_line1;

LB = -100000*ones(1,24);  % Lower bound of variables
UB = 350000*ones(1,28);
LB(25:28) = 0;

x = optimvar('x', 1, 28, 'LowerBound', LB, 'UpperBound', UB);
[total_cost, pgrid, p_pv, p_wt dccc_bat] = hazine2(x,Depot_demand,demandT1,demandT2,curi_demand_line1); 
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
prob.Constraints.cons8 = pgrid <= max(curi_demand_line1);
prob.Constraints.cons9 = p_pv - x(1:24) + p_wt <= 1.1*max(Pl);
x0.x = UB;
options = optimoptions(@fmincon,'MaxIterations',15000);
[sol,fval] = solve(prob,x0,'Options',options);
X = sol.x;
SOC(1) = X(28) / X(25);
for i=1:24
    SOC(i+1) = SOC(i) + X(i)/ X(25);
end
SOC(1) = SOC(2);
Pbat = X(1:24);
[total_cost pgrid p_pv p_wt dccc_bat] = hazine2(X,Depot_demand,demandT1,demandT2,curi_demand_line1); 
% PotherGen(PotherGen<0) = 0;

constraint4 = SOC<0.95;
D2G = pgrid.*(pgrid<0);
PV2G = D2G - 0.97*Pbat.*(Pbat<0);
G2D = pgrid.*(pgrid>0);
G2D(1) = Pl(1);
R2L = -(p_pv + p_wt - Pbat); R2L = R2L.*(R2L<0);
bar(G2D)
hold on
b = bar(R2L,'g');
hold on
bar(D2G,'r')
hold on
plot(SOC*X(25),'r','LineWidth',2)
hold on
plot(p_pv,'g-.','LineWidth',2)
hold on
% bar(PotherGen)
% plot(PotherGen,'LineWidth',2)
plot(curi_demand_line1,'LineWidth',2)
plot(curi_demand_line1 + Depot_demand + demandT1 + demandT2,'LineWidth',2)
bar(Depot_demand + demandT1 + demandT2,'y')
hold on
legend('Buy from Grid','R2L','Sell to grid','ESS','PV','City Demand','Total Load','Bus Demand')

disp(['Total Profit = ',num2str(total_cost)])
disp(['---Const 1---'])
disp(['p_pv - x(1:24) + p_wt  <= Pl: --- ','Yes'])
disp(['---Const 2---'])
disp(['p_pv - x(1:24) + p_wt <= 1.1*max(Pl): --- ','Yes'])
disp(['--------Other Output------'])
disp(['Capacity of PV: ',num2str(X(27)), ' with Max of 100000'])
disp(['Capacity of Bat: ',num2str(X(25)), ' with Max of 100000'])
disp(['Capacity of Wind: ',num2str(X(26)), ' with Max of 100000'])
% disp(['Max of Other generation: ',num2str(max(PotherGen))])
disp(['Max of City Demand: ',num2str(max(curi_demand_line1))])
fontsize = 12;
ylabel('Power (kW)')
xlabel('Time (h)')
set(gca,'fontsize',fontsize,'fontweight','b')
