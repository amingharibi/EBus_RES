function [SpeedMat,StNumMat,SlopeMat,Q_partMat,TerminalSeenMat] = Preprocess_reverse(stations,v,L)
load PassengersMatrix.mat;
station_number_go = stations(2:end,1); station_number_back = [station_number_go(end)+1:2*station_number_go(end)]';
Laststations = [station_number_go(end),station_number_back(end)];


station_dis_go = stations(2:end-1,3); distance_cumsum_go = cumsum(station_dis_go);
station_dis_back = flip(stations(2:end-1,3)); distance_cumsum_back = cumsum(station_dis_back);

stopTime = 10;
TerminalTime = 12*60;

slopes_go = stations(2:end-1,end); slopes_back = -flip(stations(2:end-1,end));
bus_number = [6,7];

T = numel(v);

s = 0;
v_new = v;

station_checked = zeros(1,20000);
timer = 2;

Q_part = 0*v_new;

passenger_stations = gradient(PassMat{1, 1}{1, 1});
passenger_stations = nonzeros(passenger_stations);
%% Back
TerminalSeen = 0;

for k=1:9
    count = 1;
    if (timer+1) + 3000 < T
        for i=timer:T
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
                        if abs(v_new(i-j) - v_new(i-j-1)) < 10
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
    %% GO
    count = 1;
    if (i3+1) + 2000 < T
        for i=i3+1:T
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
                        if abs(v_new(i-j) - v_new(i-j-1)) < 2.1
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
                        if abs(v_new(i2+j2) - v_new(i2+j2+1)) < 2
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
if L-numel(v_new) > 0
    TerminalSeen = [TerminalSeen;zeros(L-numel(v_new),1)];
    v_new = [v_new;zeros(L-numel(v_new),1)];
    station_checked = [station_checked;zeros(L-numel(station_checked),1)];
    slopes_Mat = [slopes_Mat;zeros(L-numel(slopes_Mat),1)];
    Q_part = [Q_part;zeros(L-numel(Q_part),1)];
end

v_new_extend = repmat(v_new,[1,bus_number(2)]);
station_checked_extend = repmat(station_checked',[1,bus_number(2)]);
slopes_Mat_extend = repmat(slopes_Mat',[1,bus_number(2)]);
Q_part_extend = repmat(Q_part,[1,bus_number(2)]);
TerminalSeen_extend = repmat(TerminalSeen,[1,bus_number(2)]);

timing = 10;

for i=1:bus_number(2)
    TerminalSeenMat{i} = [zeros((i-1)*timing*60,1);TerminalSeen_extend(1:size(TerminalSeen_extend,1) - (i-1)*timing*60,i)];
    SpeedMat{i} = [zeros((i-1)*timing*60,1);v_new_extend(1:size(v_new_extend,1) - (i-1)*timing*60,i)];
    StNumMat{i} = [zeros((i-1)*timing*60,1);station_checked_extend(1:size(station_checked_extend,1) - (i-1)*timing*60,i)];
    SlopeMat{i} = [zeros((i-1)*timing*60,1);slopes_Mat_extend(1:size(slopes_Mat_extend,1) - (i-1)*timing*60,i)];
    Q_partMat{i} = [zeros((i-1)*timing*60,1);Q_part_extend(1:size(Q_part_extend,1) - (i-1)*timing*60,i)];

end


end