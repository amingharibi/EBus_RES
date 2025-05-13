function [PassMat] = PassMatfiller(passenger, i_go, StNumMat,Lasstations,Bcapacity)

% passenger = passenger1;
% passenger(passenger>80) = 80;
% 
% for i=1:size(StNumMat,2)
%     PassMat{i}(1) = passenger(1);
%     for j=2:size(StNumMat{i},1)
%         PassMat{i} = [PassMat{i};PassMat{i}(j-1)];
%         if StNumMat{i}(j)~=0 && StNumMat{i}(j)~=StNumMat{i}(j-1)
%             PassMat{i}(j) = passenger(j);
%         elseif StNumMat{i}(j) == Lasstations(1) || StNumMat{i}(j) == Lasstations(2)
%             PassMat{i}(j) = 0;
%         end
%     end
% end

% for i=1:size(StNumMat,2)
%     first_index = find(StNumMat{i} ~= 0);
%     for j=first_index:size(StNumMat{i},1)
%         if StNumMat{i}(j) > Lasstations(1)
%             station = StNumMat{i}(j) - Lasstations(1);
%         else 
%             station = StNumMat{i}(j);
%         end
%         if StNumMat{i}(j) ~= StNumMat{i}(j) + 1
%             PassMat{i}(j) = PassMat2(data, station-1, j, Lasstations(1), basedTrans, Trans);
%         end
%         if StNumMat{i}(j) == Lasstations(1) || StNumMat{i}(j) == Lasstations(2)
%             PassMat{i}(j) = 0;
%         end
%     end
% end
newM2_ = zeros(1,size(StNumMat{1},1));
gain = i_go / 3600;
k = 1;
for i=1:size(StNumMat,2)
    PassMat{i}(1) = 0;
    for j=2:size(StNumMat{i},1)
        if j > 1
            newM2_(j) = newM2_(j-1);
        end
        if StNumMat{i}(j) > Lasstations(1)
            station = StNumMat{i}(j) - Lasstations(1) + 1;
        else 
            station = StNumMat{i}(j);
        end
        PassMat{i} = [PassMat{i};PassMat{i}(j-1)];
        if StNumMat{i}(j) == Lasstations(1) || StNumMat{i}(j) == Lasstations(2)
            k = j;
            PassMat{i}(j) = 0;
        elseif StNumMat{i}(j)~=0 && StNumMat{i}(j)~=StNumMat{i}(j-1)
            [PassMat{i}(j),newM] = PassMat2(passenger, gain, station-1, k, Lasstations(1), Bcapacity);
%             newM2 = max(newM,newM2);
            newM2_(j) = newM;
            
        end

    end
end
% idx = find(newM2_>0);
% newM2_(newM2_==0) = newM2_(idx(1));
% norm_newM2 = newM2_/max(newM2_);
% for i=1:numel(PassMat)
%     for j=1:numel(PassMat{1})
%         PassMat{i}(j) = norm_newM2(j)*PassMat{i}(j);
%     end
% end
% data = xlsread('پروفیل_تقاضا.xlsx',2);
% passenger = data(1:end-1,3) * (Trans/basedTrans) ;
% passenger = passenger / 62;
% passenger (1) = mean([passenger(1),passenger(2)]);
% for i=2:numel(passenger)-1
%     passenger(i) = mean([passenger(i) , passenger(i-1) , passenger(i+1)]);
% end
% passenger(end) = mean([passenger(end-1) , passenger(end)]);
% hour = floor(second/3600) + 7;
% tmp = passenger(hour)/(numstation-1)*3/2;
% a = tmp*ones((numstation-1),1);
% b = cumsum(a);
% c = round(b(1:(numstation-1)/2));
% c = [c; flipud(c)];
% tmp2 = passenger(hour)/(numstation-1)*1/2;
% c(numel(c)/2 + 1:end) = round(c(numel(c)/2 + 1:end) + tmp - tmp2);
% c(station)
end