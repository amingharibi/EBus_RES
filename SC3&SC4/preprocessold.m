clc;clear all;close all
v = xlsread('drivingCycle203.csv');
stations = xlsread('203.xlsx');
load passenger.mat;
load PassengersMatrix.mat;

passenger = passenger';
% passenger = passenger * 2;
% passenger =max(passenger)*[0,0,0,0,0,0.27,0.87,1,0.5,0.25,0.2,0.3,0.49,0.36,0.29,0.29,0.38,0.62,0.58,0.26,0.13,0.12,0.14,0.05]';
passenger = [0	0	0	0	0	196	661	886	577	441	428	445	512	519	481	517	719	885	803	435	208	164	213	99]';
station_number_go = stations(2:end,1); station_number_back = (station_number_go(end)+1:2*station_number_go(end))';
Laststations = [station_number_go(end),station_number_back(end)];

station_dis_go = stations(2:end-1,3); distance_cumsum_go = cumsum(station_dis_go);
station_dis_back = flip(stations(2:end-1,3)); distance_cumsum_back = cumsum(station_dis_back);

stopTime = 10;
TerminalTime = 12*60;
RestTim = TerminalTime;

slopes_go = stations(2:end-1,end); slopes_back = -flip(stations(2:end-1,end));
bus_number = [6,7];

passenger = passenger / sum(bus_number(1:2));
passenger (1) = mean([passenger(1),passenger(2)]);
for i=2:numel(passenger)-1
    passenger(i) = mean([passenger(i) , passenger(i-1) , passenger(i+1)]);
end
passenger(end) = mean([passenger(end-1) , passenger(end)]);


T = numel(v);


s = 0;
v_new = v;
station_checked = zeros(1,15000);
timer = 2;
Q_part = 0*[v_new;v_new];

passenger_stations = gradient(PassMat{1, 1}{1, 1});
passenger_stations = nonzeros(passenger_stations);
%% Go
TerminalSeen = 0;
for k=1:9
    count = 1;
    if (timer+1) + 3000 < T
        for i=timer:T
            s = s + v(i);
            slopes_Mat(i) = slopes_go(count);
            if s >= distance_cumsum_go(count)
                stopTime = abs(passenger_stations(count)) * 2;
                Q_part(i:i+stopTime-1) = 3.942;
                station_checked(i:i+stopTime-1) = station_number_go(count);
                count = count + 1;
                v_new(i:i+stopTime-1) = 0;
                j = 1;
                flag = 1;
                if v_new(i-1) ~= 0
                    while flag
                        v_new(i-j) = j;
                        if abs(v_new(i-j) - v_new(i-j-1)) < 2
                            flag = 0;
                        end
                        j = j+1;
                    end
                end
                i2 = i+stopTime-1;
                j2 = 1;
                flag2 = 1;
                if v_new(i2+1) ~= 0
                    while flag2
                        v_new(i2+j2) = j2;
                        if abs(v_new(i2+j2) - v_new(i2+j2+1)) < 1.9
                            flag2 = 0;
                        end
                        j2 = j2+1;
                    end
                end
            end
            if s >= max(distance_cumsum_go)
                TerminalSeen(i+1) = 1;
                TerminalTime = TerminalTime + abs(passenger_stations(count)) * 2;
                Q_part(i:i+round(normrnd(150,20))-1) = 3.942;
                station_checked(i+1:i+TerminalTime) =  station_number_go(count);
                s = 0;
                v_new(i:i+TerminalTime) = 0;
                i3 = i+TerminalTime;
                slopes_Mat(i:i+TerminalTime) = 0;
                if k==1
                    i_go = i3;
                end
                j = 1;
                flag = 1;
                if v_new(i-1) ~= 0
                    while flag
                        v_new(i-j) = j;
                        if abs(v_new(i-j) - v_new(i-j-1)) < 1.3
                            flag = 0;
                        end
                        j = j+1;
                    end
                end
                break;
            end
        end
    end
    %% Back
    count = 1;
    if (i3+1) + 3000 < T
        for i=i3+1:T
            s = s + v(i);
            slopes_Mat(i) = slopes_back(count);
            if s >= distance_cumsum_back(count)
                stopTime = abs(passenger_stations(count)) * 2;
                Q_part(i:i+stopTime-1) = 3.942;
                station_checked(i:i+stopTime-1) = station_number_back(count);
                count = count + 1;
                v_new(i:i+stopTime-1) = 0;
                j = 1;
                flag = 1;
                if v_new(i-1) ~= 0
                    while flag
                        v_new(i-j) = j;
                        if abs(v_new(i-j) - v_new(i-j-1)) < 2
                            flag = 0;
                        end
                        j = j+1;
                    end
                end
                i2 = i+stopTime-1;
                j2 = 1;
                flag2 = 1;
                if v_new(i2+1) ~= 0
                    while flag2
                        v_new(i2+j2) = j2;
                        if abs(v_new(i2+j2) - v_new(i2+j2+1)) < 1.8
                            flag2 = 0;
                        end
                        j2 = j2+1;
                    end
                end
            end
            if s >= max(distance_cumsum_back)
                TerminalSeen(i) = 1;
                TerminalTime = TerminalTime + abs(passenger_stations(count)) * 2;
                Q_part(i:i+round(normrnd(150,20))-1) = 3.942;
                s = 0;
                station_checked(i:i+TerminalTime) =  station_number_back(count+1);
                slopes_Mat(i:i+TerminalTime) = 0;
                v_new(i:i+TerminalTime) = 0;
                i3 = i+TerminalTime;
                j = 1;
                flag = 1;
                if v_new(i-1) ~= 0
                    while flag
                        v_new(i-j) = j;
                        if abs(v_new(i-j) - v_new(i-j-1)) < 1.3
                            flag = 0;
                        end
                        j = j+1;
                    end
                end
                break;
            end
        end
    end
    timer = i3 + 1;
    v_new_fake =  v_new(1:i3);
    station_checked = station_checked(1:i3);
    Q_part_ = Q_part(1:i3);
    Q_part_cycles{k} = Q_part;
    driving_cycles{k} = v_new_fake;
    station_checked_cycles{k} = station_checked;
    TerminalSeen_cycles{k} = TerminalSeen;
end
TerminalSeen = TerminalSeen_cycles{end}';
Q_part = Q_part_cycles{end};
v_new = driving_cycles{end};
station_checked = station_checked_cycles{end};
TerminalSeen(end+1:end+1+numel(station_checked)-numel(TerminalSeen))=0;

v_new_extend = repmat(v_new,[1,bus_number(1)]);
station_checked_extend = repmat(station_checked',[1,bus_number(1)]);
slopes_Mat_extend = repmat(slopes_Mat',[1,bus_number(1)]);
Q_part_extend = repmat(Q_part,[1,bus_number(1)]);
TerminalSeen_extend = repmat(TerminalSeen,[1,bus_number(1)]);

timing = 10;

for i=1:bus_number(1)
    TerminalSeenMat{i} = [zeros((i-1)*timing*60,1);TerminalSeen_extend(1:size(TerminalSeen_extend,1) - (i-1)*timing*60,i)];
    SpeedMat{i} = [zeros((i-1)*timing*60,1);v_new_extend(1:size(v_new_extend,1) - (i-1)*timing*60,i)];
    StNumMat{i} = [zeros((i-1)*timing*60,1);station_checked_extend(1:size(station_checked_extend,1) - (i-1)*timing*60,i)];
    SlopeMat{i} = [zeros((i-1)*timing*60,1);slopes_Mat_extend(1:size(slopes_Mat_extend,1) - (i-1)*timing*60,i)];
    Q_partMat{i} = [zeros((i-1)*timing*60,1);Q_part_extend(1:size(Q_part_extend,1) - (i-1)*timing*60,i)];
end

Bcapacity = 270;

PassMat = PassMatfiller(passenger, i_go, StNumMat, Laststations, Bcapacity);

L = size(v_new_extend,1);

[SpeedMat_reverse,StNumMat_reverse,SlopeMat_reverse,Q_partMat_reverse,TerminalSeenMat_reverse] = Preprocess_reverse(stations,v,L);

PassMat_reverse = PassMatfiller(passenger, i_go, StNumMat_reverse, Laststations, Bcapacity);

for ii=1:size(StNumMat,2)
    Ix = find((StNumMat{ii}~=0));
    StNumMat{ii}(Ix(end)) = -1;
end

for ii=1:size(StNumMat_reverse,2)
    Ix = find((StNumMat_reverse{ii}~=0));
    StNumMat_reverse{ii}(Ix(end)) = -1;
end

Q_partMat_overall{1} = Q_partMat;
Q_partMat_overall{2} = Q_partMat_reverse;


SpeedMat_overall{1} = SpeedMat;
SpeedMat_overall{2} = SpeedMat_reverse;

StNumMat_overall{1} = StNumMat;
StNumMat_overall{2} = StNumMat_reverse;

SlopeMat_overall{1} = SlopeMat;
SlopeMat_overall{2} = SlopeMat_reverse;

PassMat_overall{1} = PassMat;
PassMat_overall{2} = PassMat_reverse;

TerminalSeenMat_overall{1} = TerminalSeenMat;
TerminalSeenMat_overall{2} = TerminalSeenMat_reverse;

TerminalSeenMat = TerminalSeenMat_overall;
SpeedMat = SpeedMat_overall;
StNumMat = StNumMat_overall;
SlopeMat = SlopeMat_overall;
PassMat = PassMat_overall;
QpartMat = Q_partMat_overall;




for i=1:floor(L/3600)
    min_bus_1(i) = ceil(max(PassMat{1}{1}((i-1)*3600+1:i*3600)));
    min_bus_2(i) = ceil(max(PassMat{1}{1}((i-1)*3600+1:i*3600)));
    %     min_bus_1(i) = ceil(passenger(i))/2;
    %     min_bus_2(i) = ceil(passenger(i))/2;
end

minBusN(:,1) = 6 + 0*ceil(min_bus_1*bus_number(1)/Bcapacity);
minBusN(:,2) = 6 + 0*ceil(min_bus_2*bus_number(2)/Bcapacity);
% minBusN(:,1) = 6;
% minBusN(:,2) = 7;
% minBusN_timing = round((i_go-TerminalTime)/(timing*60));
% minBusN(minBusN(:,1)<minBusN_timing,1) = minBusN_timing;
% minBusN(minBusN(:,2)<minBusN_timing,2) = minBusN_timing;

Time = 1:L;

for i=1:size(StNumMat,2)
    for j=1:size(StNumMat{i},2)
        numTermialSeen{i}{j} = numel(find(TerminalSeenMat{i}{j}==1));
        L_tmp = size(StNumMat{i}{j},1);
        FlagMat{i}{j} = ones(L_tmp,1);
        ChrgFlag{i}{j} = zeros(L_tmp,1);
    end
end

Bus_Num_Mat = [bus_number(1);bus_number(2)];

tariff = 0;
tariff(1:12*3600) = 0.6;
tariff(12*3600:13*3600) = 1;
tariff(13*3600:16*3600) = 1.6;
tariff(16*3600:17*3600) = 1;
% tariff = tariff * 0 + 0.45;
tariff = tariff*0.2;

Temp_hourly = [59,61,64,70,75,77,82,84,86,88,88,88,81,75,72,70,68];
Temp_hourly = (Temp_hourly-32) * (5/9);
amb_temp = [];
for i=1:numel(Temp_hourly)-1
    amb_temp(3600*(i-1)+1:3600*i) = linspace(Temp_hourly(i),Temp_hourly(i+1),3600);
end

save integrated_data2 SpeedMat StNumMat SlopeMat PassMat L minBusN Time FlagMat ChrgFlag Bus_Num_Mat Laststations RestTim tariff QpartMat amb_temp numTermialSeen TerminalSeenMat

