function J = Cost(inputs)

x_ = round(inputs);
offset = 21600;
load integrated_data2.mat
minBusN = [5,6;5,6;6,7;5,5;4,5;4,4;4,5;4,5;4,5;4,5;5,6;6,7;6,7;6,6];

t_lb = 22;
t_s = 25;
t_ub = 28;
UpperSOC = 95;
LimitSOC = 5;
charger_cost = [62000,87500,175000,247500];
Charg_rates = [124,180,350,450];

DepotCharger_cost = [9000 24500 41000 62000];
DepotCharger_rate = [20 40 60 120];

charge_choose(1) = x_(4);
charge_choose(2) = x_(5);

deposit_energy = 78; deposit_SOC = 100*deposit_energy/x_(1);
X = x_(6:end);
counter = 0;
for i=1:size(StNumMat,2)
    for j=1:size(StNumMat{i},2)
        numberofTerminalseen{i}{j} = 1;
        Charge{i}{j} = X(counter+1:numTermialSeen{i}{j}+counter);
        charge_cumsum = 0;
        for z = 1:numTermialSeen{i}{j}
            if rem(z,2) == 0
                charge_cumsum = charge_cumsum + sum(Charge{i}{j}*Charg_rates(x_(5)));
            else
                charge_cumsum = charge_cumsum + sum(Charge{i}{j}*Charg_rates(x_(4)));
            end
        end
        if (2*deposit_energy + 750 - charge_cumsum)/7 > DepotCharger_rate(x_(end))
            J=inf;
            return;
        end
        
        counter = counter + numTermialSeen{i}{j};
        numberofTerminalseen{i}{j} = 0;
        Cost_charge{i}{j} = 0;
        ChargeEnergy{i}{j} = 0;
    end
end

%% Load data
% numTermialSeen TerminalSeenMat
%% Parameters
mv = 40e3*1.25; % vehicle curb mass (instead of 50 ton for long term mass while Battery pack weight will be reduced)
Cd = 0.6; % drag coef.
Af = 6.1; % vehicle eff. area
rho = 1.2258; % air density
g = 9.81; % gravity
effMot = 0.85; % total drivetrain efficiency (Motor)
effReg = 0.67; % total drivetrain efficiency (Regen)
f = 0.01;
auxPower = 6000;
fancoil = 5724;
%% Battery pack

nomEnergy = x_(1);
PantoCstTime = 0; % 15+15+15 % Equal Time to 1) arrive bus in to Panto charger at proposed station 2) install 3) uninstall Panto 4) Battery cooling power dissipation equal time (sec)
ChrgLim23 = x_(2);
ChrgLim46 = x_(3);

mv = mv + 14.29*(nomEnergy-792); % vehicle curb mass + extra battery pack (160 Wh/kg EVE)

Charger_rate(1) = Charg_rates(charge_choose(1));
Charger_rate(2) = Charg_rates(charge_choose(2));

chGain = nomEnergy./Charger_rate;
chGain(chGain<1) = 1;

status{1} = ones(1,Bus_Num_Mat(1));
status{2} = ones(1,Bus_Num_Mat(2));
minBusN = [minBusN1' minBusN2'];
minBusN(end+1,:) = minBusN(end,:) - 1;
E = nomEnergy;
discount = -0.04613;
discount2 = 3.28/100;
BatteyCost = 188;
CC = BatteyCost*nomEnergy;
T_bus = 25;

%% Main Loop
i = 0;
ULSOC = UpperSOC; % percent
ULTime(1) = interp1([0,100],[0,3600*chGain(1)],ULSOC,'linear','extrap');
ULTime(2) = interp1([0,100],[0,3600*chGain(2)],ULSOC,'linear','extrap');

ChdargePerCycle = Charg_rates(charge_choose).*(RestTim/3600); 

while i<L-1
    i = i + 1;
    if rem(i,3600) == 0
        i;
    end
    for j = 1:numel(Bus_Num_Mat)
        for k = 1:Bus_Num_Mat(j)
            if (StNumMat{j}{k}(i)==2 || StNumMat{j}{k}(i)==33)
                if Energy{j}{k}(i-1)<30
                    J=inf;
                    return;
                end
            end
            if i==1
                gamma(i) = -1;
                T_bus(2) = 25;
                Energy{j}{k}(i) = nomEnergy; % battery energy (Wh)
                SOC{j}{k}(i) = interp1([0,nomEnergy],[0,100],Energy{j}{k}(i),'linear','extrap');
            
            elseif ~ChrgFlag{j}{k}(i)
                Vb = SpeedMat{j}{k}(i); % vehicle speed (m/s)
                dVb = (SpeedMat{j}{k}(i)-SpeedMat{j}{k}(i-1));
                slope = SlopeMat{j}{k}(i); % slope (%)
                theta = atan(slope/100); % slope (rad)
                mp = 76*PassMat{j}{k}(i); % passenger mass (kg)
                m = (1.02*mv+mp); % vehicle mass + passenger mass (kg)
                Fg = m*g*sin(theta); % gravity force
                Fr = m*g*f*cos(theta); % rolling resistance force
                Fa = 0.5*Cd*rho*Af*Vb^2; % aerodynamic force
                if PassMat{j}{k}(i)>50
                    dVb(dVb>1) = 1;
                end
                Fi = m*dVb; % inertia force
                Ftb = Fi+Fa+Fr+Fg; % traction/brake force

                delta_temp = (T_bus(i) - amb_temp(i));
                Q_t_1{j}{k}(i) = delta_temp / (0.0051 + 1/(271.74*(10.45-SpeedMat{j}{k}(i) + 10*sqrt(SpeedMat{j}{k}(i)))));
                Q_t_2{j}{k}(i) = 0.556*delta_temp;
                Q_t_3{j}{k}(i) = QpartMat{j}{k}(i);
                Q_t_4{j}{k}(i) = 1.8*PassMat{j}{k}(i);
                t_bus = T_bus(i);

                if i>1
                    if (t_bus < t_lb) || (gamma(i-1) == 1 && t_bus < t_s)
                        gamma(i) = +1;
                    elseif (t_bus > t_ub) || (gamma(i-1) == -1 && t_bus > t_s)
                        gamma(i) = -1;
                    else
                        gamma(i) = 0;
                    end
                end
                Q_HVAC{j}{k}(i) = gamma(i) * 60000;
                T_bus(i+1) =  + ( ((Q_HVAC{j}{k}(i) - Q_t_1{j}{k}(i) - Q_t_2{j}{k}(i) - Q_t_3{j}{k}(i) + Q_t_4{j}{k}(i))/(1000*131.5*1.005)) + t_bus);

                E_t = abs(Q_HVAC{j}{k}(i)) / (1.9000);


                if (StNumMat{j}{k}(i)==Laststations(1)) || (StNumMat{j}{k}(i)==Laststations(2))
                    auxPower = 0;
                end
                Pax = E_t + (auxPower+fancoil);
                Pax_{j}{k}(i) = Pax/1e3/3600;
                switch sign(Ftb*Vb)
                    case -1
                        ConsPower{j}{k}(i) = FlagMat{j}{k}(i)*(Vb*Ftb*effReg+Pax)/1e3; % Consumed energy in 1 second (kWh)
                        ConsEnergy{j}{k}(i) = FlagMat{j}{k}(i)*(Vb*Ftb*effReg+Pax)/1e3/3600; % Consumed energy in 1 second (kWh)
                    otherwise
                        ConsPower{j}{k}(i) = FlagMat{j}{k}(i)*(Vb*Ftb/effMot+Pax)/1e3; % Consumed energy in 1 second (kWh)
                        ConsEnergy{j}{k}(i) = FlagMat{j}{k}(i)*(Vb*Ftb/effMot+Pax)/1e3/3600; % Consumed energy in 1 second (kWh)
                end
                Energy{j}{k}(i) = Energy{j}{k}(i-1)-ConsEnergy{j}{k}(i); % battery energy (kWh)
                SOC{j}{k}(i) = interp1([0,nomEnergy],[0,100],Energy{j}{k}(i),'linear','extrap');

                if SOC{j}{k}(i)<LimitSOC
                    J = inf; return;
                end

            end
            if TerminalSeenMat{j}{k}(i) == 1
                numberofTerminalseen{j}{k} = numberofTerminalseen{j}{k} + 1;
            end
            flag_terminal = 0;
            if StNumMat{j}{k}(i)==Laststations(1)
                flag_terminal = 1;
            elseif StNumMat{j}{k}(i)==Laststations(2)
                flag_terminal = 2;
            end
            if (StNumMat{j}{k}(i)==Laststations(1) || StNumMat{j}{k}(i)==Laststations(2)) && ...
                    Charge{j}{k}(numberofTerminalseen{j}{k}) && ~ChrgFlag{j}{k}(i) && ~ChrgFlag{j}{k}(i-RestTim) && ...
                    (ULSOC - SOC{j}{k}(i))*nomEnergy/100 >= ChdargePerCycle(flag_terminal)

                c = 0;
                for m1 = 1:numel(Bus_Num_Mat)
                    for m2 = 1:Bus_Num_Mat(m1)
                        c = c+1;
                        if i + 1 <= numel(ChrgFlag{m1}{m2})
                            chflag = ChrgFlag{m1}{m2}(i+1);
                            iCh(c) = chflag;
                        else
                            chflag = ChrgFlag{m1}{m2}(end);
                            iCh(c) = chflag;
                        end
                    end
                end

                SumChrg23 = sum(iCh==Laststations(1)); SumChrg46 = sum(iCh==Laststations(2));
                SumChrg23_(i) = SumChrg23; SumChrg46_(i) = SumChrg46;
                if (StNumMat{j}{k}(i)==Laststations(1) && SumChrg23<ChrgLim23) || (StNumMat{j}{k}(i)==Laststations(2) && SumChrg46<ChrgLim46)                     
                    s = 0;
                    for i1 = 1:Bus_Num_Mat(j)
                        if i + 1 <= numel(ChrgFlag{j}{i1})
                            s = s+~(~ChrgFlag{j}{i1}(i+1));
                        else
                            s = s+~(~ChrgFlag{j}{i1}(end));
                        end
                    end
                    FlagMinBus = 0;
                    if round(i/3600) <= 16
                        if Bus_Num_Mat(j)-s>=minBusN(ceil(i/3600),j); FlagMinBus = 1; end
                    else
                        if Bus_Num_Mat(j)-s>=minBusN(end,j); FlagMinBus = 1; end
                    end
                    if StNumMat{j}{k}(i)==Laststations(1)
                        aa = find(StNumMat{j}{k}~=Laststations(1));
                        aa = aa(aa>i);
                        allowChtime = aa(1)-i;
                    else
                        aa = find(StNumMat{j}{k}~=Laststations(2));
                        aa = aa(aa>i);
                        allowChtime = aa(1)-i;
                    end
                    allowChtime = min(allowChtime,RestTim);%2*RestTim - sum(ChrgFlag{j}{k}(:)>0);
                    if ~FlagMinBus && SOC{j}{k}(i)< LimitSOC
                        J = inf;
                        return
                    elseif  FlagMinBus && (allowChtime>0)
                        if StNumMat{j}{k}(i)==Laststations(1)
                            ULTime_tmp = ULTime(1); chGain_tmp = chGain(1);
                        else
                            ULTime_tmp = ULTime(2); chGain_tmp = chGain(2);
                        end
                        PantoChTime = PantoCstTime+ceil(ULTime_tmp-interp1([0,100],[0,3600*chGain_tmp],SOC{j}{k}(i),'linear','extrap'));
                        PantoChTime(PantoChTime>RestTim) = RestTim;
                        PantoChTime = min(PantoChTime,allowChtime);
                        ULEnergy = interp1([0,60*chGain_tmp],[0,nomEnergy],PantoChTime/60,'linear','extrap');
                        ChargeEnergy{j}{k}(i+1:i+PantoChTime) = (ULEnergy./(PantoChTime));
                        Energy{j}{k}(i+1:i+PantoChTime) = interp1([i,i+PantoChTime],[Energy{j}{k}(i),Energy{j}{k}(i) + ULEnergy],i+1:i+PantoChTime);
                        Cost_charge{j}{k}(i+1:i+PantoChTime) = (ULEnergy./(PantoChTime)) * tariff(i+1:i+PantoChTime);
                        SOC{j}{k}(i+1:i+PantoChTime) = interp1([0,nomEnergy],[0,100],Energy{j}{k}(i+1:i+PantoChTime),'linear','extrap');
                        ChrgFlag{j}{k}(i+1:i+PantoChTime) = StNumMat{j}{k}(i)*ones(PantoChTime,1);
                        FlagMat{j}{k}(i+1:i+PantoChTime) = 0;
                    end
                end
            end
        end
    end
end

BatCost = [];
depotcharge = 0;
charge_cost = [];
Charge_Energy_ = [];
depotcharge_ = [];
SumChrg23_ = 0;
SumChrg46_ = 0;
Demand23_ = 0;
Demand46_ = 0;

maxL = 0;
for i=1:size(Cost_charge,2)
    for j=1:size(Cost_charge{i},2)
        if size(ChargeEnergy{i}{j},2) >= maxL
            maxL = size(ChargeEnergy{i}{j},2);
        end
    end
end
ChargeEnergy_tmp = ChargeEnergy;
for i=1:size(Cost_charge,2)
    for j=1:size(Cost_charge{i},2)
        if SOC{i}{j}(end) < deposit_SOC + LimitSOC
            J = inf;
            % return
        end
        ChargeEnergy_tmp{i}{j} = [ChargeEnergy_tmp{i}{j},zeros(1,maxL-size(ChargeEnergy{i}{j},2))];
        SumChrg23_ = SumChrg23_ + sign(ChrgFlag{i}{j}(1:L-1).*(StNumMat{i}{j}(1:L-1)==Laststations(1)));
        SumChrg46_ = SumChrg46_ + sign(ChrgFlag{i}{j}(1:L-1).*(StNumMat{i}{j}(1:L-1)==Laststations(2)));
        
        Demand23_ = Demand23_ + ChargeEnergy_tmp{i}{j}(1:maxL).*(StNumMat{i}{j}(1:maxL)==Laststations(1))';
        Demand46_ = Demand46_ + ChargeEnergy_tmp{i}{j}(1:maxL).*(StNumMat{i}{j}(1:maxL)==Laststations(2))';

        charge_cost = [charge_cost, sum(Cost_charge{i}{j})];
        Charge_Energy = sum(ChargeEnergy{i}{j}) + nomEnergy * (UpperSOC-SOC{i}{j}(end))/100;
        Charge_Energy_ = [Charge_Energy_ sum(ChargeEnergy{i}{j})];
        depotcharge = depotcharge + nomEnergy * (UpperSOC-SOC{i}{j}(end))/100;
        depotcharge_ = [depotcharge_ nomEnergy * (UpperSOC-SOC{i}{j}(end))/100];
        BatCost = [BatCost deltaE(nomEnergy,Charge_Energy,discount,CC)];
    end
end
for i=1:13
    if depotcharge_(i)/7.5>DepotCharger_rate(x_(end))
        J=inf;
        return
    end
end
Charge_Energy = sum(Charge_Energy);
BatCost = sum(BatCost);
Energy_cost = sum(charge_cost);
ChargerCostRF = chargerRF(discount2,charger_cost(charge_choose));
ChargerCostRF_depot = chargerRF(discount2,DepotCharger_cost(x_(end)));

% mean((depotcharge_ + Charge_Energy_)./nomEnergy)

Depot_demand = zeros(1,24); Depot_demand(1,1:5) = sum(depotcharge_ / 7);
Depot_demand(1,22:24) = sum(depotcharge_ / 7);

demandT2 = zeros(1,24); demandT1 = zeros(1,24); 

SumChrg46__ = [SumChrg46_;zeros(3600*15 - size(SumChrg46_,1),1)];
SumChrg23__ = [SumChrg23_;zeros(3600*15 - size(SumChrg23_,1),1)];

data_hours46 = reshape(Charger_rate(2) * SumChrg46__, 3600, [])'; 
data_hours23 = reshape(Charger_rate(1) * SumChrg23__, 3600, [])'; 

hourly_mean46 = mean(data_hours46, 2);
hourly_mean23 = mean(data_hours23, 2);

demandT1(1,6:20) = hourly_mean23;
demandT2(1,6:20) = hourly_mean46;
%% Demands
T1DemandCost = max(Charger_rate(1) * SumChrg23__) * 11.5 * 1.08; % 1.08: euro to USD
T2DemandCost = max(Charger_rate(2) * SumChrg46__) * 11.5 * 1.08;
DepotDemandCost = max(Depot_demand) * 11.5 * 1.08;

DemandCost = T1DemandCost + T2DemandCost + DepotDemandCost;
inflation = 0.0328/365;
demandcost = ((inflation*DemandCost)/ ( ((1+inflation)^30) - 1));

J = sum(ChargerCostRF.*(x_(2:3))) + Energy_cost + BatCost + depotcharge*min(tariff) + demandcost + 13*ChargerCostRF_depot;

%%
% offset = 21600;
% fontsize = 12;
% for i1=1:Bus_Num_Mat(2)
%     subplot(2,1,1)
%     plot((offset+1:(L-1)+offset)/3600,SOC{1, 2}{1, i1}(1:L-1),'LineWidth',1.5)
%     hold on
%     ylabel('SOC')
%     xlabel('Time (h)')
% end
% xlim([0,24])
% ylim([0,150])
% xline(6,'-',{'Start'},'LabelHorizontalAlignment','center','LineWidth',1.5);
% xline(((L-1)+offset)/3600,'-',{'End'},'LabelHorizontalAlignment','center','LineWidth',1.5);
% yline(30,'-',{'SOC=30'},'LabelVerticalAlignment','bottom','LineWidth',1.5);
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% offset = 21600;
% for i2=1:Bus_Num_Mat(1)
%     subplot(2,1,2)
%     plot((offset+1:(L-1)+offset)/3600,SOC{1, 1}{1, i2}(1:L-1),'LineWidth',1.5)
%     hold on
%     ylabel('SOC')
%     xlabel('Time (h)')
% end
% xlim([0,24])
% ylim([0,150])
% xline(6,'-',{'Start'},'LabelHorizontalAlignment','center','LineWidth',1.5);
% xline(((L-1)+offset)/3600,'-',{'End'},'LabelHorizontalAlignment','center','LineWidth',1.5);
% yline(30,'-',{'SOC=30'},'LabelVerticalAlignment','bottom','LineWidth',1.5);
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% % %% 
% offset = 21600;
% fontsize = 12;
% fig = figure;
% 
% X_axis = (offset+1:(L-1)+offset)/3600;
% subplot(2,1,1)
% plot(X_axis,Charger_rate(1)*SumChrg23_,'LineWidth',2)
% hold on
% ylabel('Load Demand in Terminal 1 (kW)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% subplot(2,1,2)
% plot(X_axis,Charger_rate(2)*SumChrg46_,'LineWidth',2)
% ylabel('Load Demand in Terminal 2 (kW)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% han=axes(fig,'visible','off');
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% xlabel(han,'Time (h)');
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% 
% % for article plot
% secondly_matrix_repeated = repelem(curi_demand_line1, 1, 3600);
% secondly_matrix_repeated(21600:53053+21600-1)'+Charger_rate(1)*SumChrg23_+Charger_rate(2)*SumChrg46_;
% 
% offset = 21600;
% fontsize = 12;
% fig = figure;
% 
% X_axis = (offset+1:(L-1)+offset)/3600;
% subplot(2,1,1)
% plot(X_axis,Charger_rate(1)*SumChrg23_,'LineWidth',2)
% hold on
% ylabel('Load Demand in Terminal 1 (kW)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% subplot(2,1,2)
% plot(X_axis,Charger_rate(2)*SumChrg46_,'LineWidth',2)
% ylabel('Load Demand in Terminal 2 (kW)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% han=axes(fig,'visible','off');
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% xlabel(han,'Time (h)');
% set(gca,'fontsize',fontsize,'fontweight','b')

% % 
% % Depot_demand = zeros(1,24); Depot_demand(1,1:5) = sum(depotcharge_ / 7);
% % Depot_demand(1,22:24) = sum(depotcharge_ / 7);
% % 
% % demandT2 = zeros(1,24); demandT1 = zeros(1,24); 
% % 
% % SumChrg46__ = [SumChrg46_;zeros(3600*15 - size(SumChrg46_,1),1)];
% % SumChrg23__ = [SumChrg23_;zeros(3600*15 - size(SumChrg23_,1),1)];
% % 
% % data_hours46 = reshape(Charger_rate(2) * SumChrg46__, 3600, [])'; 
% % data_hours23 = reshape(Charger_rate(1) * SumChrg23__, 3600, [])'; 
% % 
% % hourly_mean46 = mean(data_hours46, 2);
% % hourly_mean23 = mean(data_hours23, 2);
% % 
% % demandT1(1,6:20) = hourly_mean23;
% % demandT2(1,6:20) = hourly_mean46;
% % 
% Pl = demandT1 + demandT2 + Depot_demand;
% plot(Pl,'LineWidth',3)
% ylabel('Load Demand in Terminal and Depot (kW)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% xlabel('Time (h)');
% set(gca,'fontsize',fontsize,'fontweight','b')

% hold on
% plot(demandT2,'LineWidth',2)
% hold on
% plot(Depot_demand,'LineWidth',2)



% 
% %%
% tariff = 0;
% tariff(1:18*3600) = 0.6;
% tariff(18*3600:20*3600) = 1;
% tariff(20*3600:22*3600) = 1.6;
% tariff(22*3600:24*3600) = 1;
% tariff(24*3600:25*3600) = 0.6;
% tariff = tariff*0.2;
% load passenger.mat;
% 
% minBusN = [0 0;0 0; 0 0; 0 0; 0 0; 0 0; minBusN; minBusN(end,1)-2 minBusN(end,2)-2; 0 0; 0 0];
% MinBusN = zeros(86400,2);
% for i=1:24
%     Pass((i-1)*3600+1:i*3600) = passenger(i);
%     for j=1:2
%         MinBusN((i-1)*3600+1:i*3600,j) = minBusN(i,j);
%     end
% end
% 
% offset = 21600;
% fontsize = 10;
% fig = figure;
% 
% X_axis = (offset+1:(L-1)+offset)/3600;
% subplot(3,1,1)
% yyaxis left
% plot(X_axis,Charger_rate(1)*SumChrg23_,'LineWidth',2)
% ylim([0 200])
% ylabel('Demand T1 (kW)')
% hold on
% yyaxis right
% plot([1:numel(tariff)]/3600,tariff,'LineWidth',2)
% ylabel('Tariff ($)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% xlim([6,22])
% 
% subplot(3,1,2)
% yyaxis left
% plot(X_axis,Charger_rate(2)*SumChrg46_,'LineWidth',2)
% ylim([0 200])
% ylabel('Demand T2 (kW)')
% yyaxis right
% plot([1:numel(tariff)]/3600,tariff,'LineWidth',2)
% ylabel('Tariff ($)')
% set(gca,'fontsize',fontsize,'fontweight','b')
% xlim([6,22])
% 
% 
% 
% X_axis = (offset+1:(L-1)+offset)/3600;
% subplot(3,1,3)
% yyaxis left
% plot([1:numel(Pass)]/3600,Pass,'LineWidth',2)
% ylim([0 1200])
% ylabel('Passenger')
% hold on
% yyaxis right
% plot([1:numel(MinBusN(:,1))]/3600,MinBusN(:,1),'LineWidth',2,'Color','r')
% hold on 
% plot([1:numel(MinBusN(:,2))]/3600,MinBusN(:,2),'LineWidth',2,'Color','r','LineStyle','--')
% ylim([0 8])
% lgd = legend('Passenger','Forward route','Back route','Location','best');
% lgd.FontSize = 8;
% ylabel('MinBus')
% set(gca,'fontsize',fontsize,'fontweight','b')
% xlim([6,22])
% 
% han=axes(fig,'visible','off');
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% xlabel(han,'Time (h)');
% 
% set(gca,'fontsize',fontsize,'fontweight','b')
% 
% 
% end