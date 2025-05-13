function [pass,exitM] = PassMat2(passenger, gain, station, second, numstation, Bcapacity)

hour = floor(second/3600) + 7;
M = passenger(hour);%/(2*(numstation-1));
P = Bcapacity; % Capacity of E-bus
exitM = passenger(hour)*gain;
newM = gain*M; % Consider portion of cycle in 1 hour
S = numstation; % Number of stations in 1 direction

for k = 1:length(newM)
    %% Initial add passenger to just meet P
    y = zeros(S,S);
    x = zeros(1,S);
    xTot = zeros(1,S);
    yTot = zeros(1,S);
    for i = 1:S-1
        x(i) = newM/(S-1); % Passengers go to
        for j = i+1:S
            y(j,i) = x(i)/(S-i); % Passengers come out
        end
        yTot(i) = sum(y(i,:));
    end
    xTot = cumsum(x);
    yTott = cumsum(yTot);
    
    totPass = cumsum(x-yTot); totPass(end) = 0;
    normalPass = totPass*P/max(totPass);
    finalPass = round(normalPass*M/max(passenger));
end

pass = finalPass(station);
% %% Initial add passenger to just meet Bcapacity
% N = min(newM,Bcapacity); % Peak passengers of E-bus during scenario
% V = cumsum(Bcapacity/(numstation/2-1)*ones(1,floor(numstation/2-1)))';
% V = [V; flipud(V); 0];
% 
% %% Second add passenger to just meet M
% if newM>Bcapacity
%     newP = newM-Bcapacity;
%     V([2:floor(numstation/2)-2,floor(numstation/2)+2:numstation-2]) = V([2:floor(numstation/2)-2,floor(numstation/2)+2:numstation-2])+newP/length([2:floor(numstation/2)-2,floor(numstation/2)+2:numstation-2]);
% end
% 
% %% Substract passenger just to simulate passenger down
% sigma = M;
% V(5:numstation-5) = V(5:numstation-5)-abs(sigma*randn(length((5:numstation-5)),1));
% pass = V(station);

% S = 37; % Number of stations in 1 direction
% P = 80; % Capacity of E-bus
% Tp = 110; % Time takes for passing 1 direction
% T = 60; % Reference time
% newM = Tp/T*M(k); % Consider portion of cycle in 1 hour
% 
% %% Initial add passenger to just meet P
% y = zeros(S,S);
% x = zeros(1,S);
% xTot = zeros(1,S);
% yTot = zeros(1,S);
% for i = 1:S-1
%     x(i) = newM/(S-1); % Passengers go to
%     for j = i+1:S
%         y(j,i) = x(i)/(S-i); % Passengers come out
%     end
%     yTot(i) = sum(y(i,:));
% end
% xTot = cumsum(x);
% yTott = cumsum(yTot);
% 
% totPass = cumsum(x-yTot); totPass(end) = 0;
% normalPass = totPass*P/max(totPass);
% finalPass = normalPass*M(k)/maxM;
end