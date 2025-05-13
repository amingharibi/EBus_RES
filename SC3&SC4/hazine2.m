function [z, pgrid, p_pv, PotherGen, p_wt,dccc_bat NECL TMP]=hazine2(x,Pl1,Pl2,Pl3,curi_demand_line1)


inf=0.0328/365;
batprice = 139;
SVc=1.3*(1/0.453)*14.26;

f_pv=0.9;
Kt=-3.7*10e-3;
Tref=25;
v_citin=3;
v_cutout=25;
v_rated=13;
z_w=3;

wt_cc=690.41;     %$/kW
wt_sv=0;      %$/kW
wt_mc=21.45;        %$/kW
wt_l=25*365;

pv_cc=460.27;     %$/kW
pv_sv=0;      %$/kW
pv_mc=1.73;        %$/kW
pv_l=25*365;

Mu=[0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 1 1.6 1.6 1.6 1 0.6 0.6]*0.2;
% Mu = [0.380146705005575	0.389047591103809	0.393648260078634	0.394512059151458	0.390756411008744	0.385536060090370	0.373424094830116	0.353594272636582	0.340374391174227	0.330741153688164	0.327924417581128	0.334290241183029	0.328731881931812	0.323905874068423	0.323830761105569	0.326215597676193	0.330346810633179	0.320788686109970	0.320000000000000	0.320037556481427	0.322253388885629	0.331417170353853	0.347660348571093	0.366814154098938];
Pl=Pl1+Pl2+Pl3+curi_demand_line1;
P_bus=Pl1+Pl2+Pl3;

%% input renewable
windspeed=[6,6,6,6,6,6,5,5,5,5,5,5,6,6,6,5,5,6,6,6,6,6,6,6];   %bus in 150m
windspeed = windspeed * (200/100)^0.14;
solarir=[0,0,0,0,0,0,0,1,27,89,144,175,198,214,221,218,202,170,122,57,10,0,0,0]; %bus
Tamb=[16,16,15,15,15,15,14,14,14,15,17,19,21,22,23,23,23,23,22,20,19,18,17,16];
%%


c_wt=x(26); % min = 0 max = 1.1
c_pv=x(27); % min = 0 
  


for i=1:24
   
   if windspeed(i) < v_citin || windspeed(i) > v_cutout
       p_wt(i)=0;
   end
   if windspeed(i)>= v_citin && windspeed(i) <= v_rated
       p_wt(i)=c_wt*(((windspeed(i)-v_citin)/(v_rated-v_citin))^z_w);
   end      
   if windspeed(i) > v_rated && windspeed(i) <= v_cutout 
       p_wt(i)=c_wt+((0.72*c_wt-c_wt)/(v_cutout -v_rated))*(windspeed(i)-v_rated);
   end    

   p_pv(i) = c_pv*f_pv*(solarir(i)/1000)*(1+Kt*((Tamb(i)+0.0256*solarir(i))-Tref));
     
end       
pbat = 0.927*x(1:24);
pgrid = pbat - p_wt - p_pv;
z = sum(pgrid .* Mu);
PotherGen = Pl + pgrid;

 
% NECL = 0;
% dccc_bat = ((inf*((batprice*batcap*((1+inf).^NECL)) - SVc*batcap))/ (((1+inf).^NECL) - 1)); 
% 
% dcc_wt = ((inf*((wt_cc*c_wt*((1+inf)^wt_l)) - wt_sv*c_wt))/ ( ((1+inf)^wt_l) - 1));
% 
% dcc_pv = ((inf*((pv_cc*c_pv*((1+inf)^pv_l)) - pv_sv*c_pv))/ ( ((1+inf)^pv_l) - 1));
% 
% doc_wt=((inf*wt_mc*c_wt)/ ( ((1+inf)^365) - 1));
% 
% doc_pv=((inf*pv_mc*c_pv)/ ( ((1+inf)^365) - 1));

% z=2;

% J=sum(Cu)+dcc_wt+dcc_pv;
batcap = x(25);
tmp = sqrt(x(1:24).^2+0.00001);
% NECL = 215.2068 +  19.8454*sum(tmp) + -3.1941*batcap + -0.0652*batcap*sum(tmp) + 0.0105*batcap^2;
% y = sum(tmp)/1000;

% NECL = 622.1688 + 190.6727*batcap + -88.9671*y + 0.0025*batcap^2 + -9.5617*batcap*y + ...
%     4.4606*y^2 + -0.0001*batcap^2*y + 0.2021*batcap*y^2 + -0.0942*y^3 + -0.0015*batcap*y^3 + 0.0007*y^4;
NECL = (0.2*batcap)/(1*sum(tmp)*0.00015);
% dccc_bat = ((inf*((batprice*batcap*((1+inf)^NECL)) - SVc*batcap))/ (((1+inf)^NECL) - 1)); 
dccc_bat = (1*batprice/139) * 1.06*batcap * 98.0817 / NECL + 0.0081;
% NECL_ = NECL/100;
% dccc_bat = NECL_;
% NECL_ = 1;
% NECL_ = sum(tmp);
% dccc_bat = batcap * (2.6826 * NECL_^4 + -7.8399 * NECL_^3 + 8.6803 * NECL_^2 + -4.5134 * NECL_^1 + 1.1012);
TMP = sum(tmp);
dcc_wt = ((inf*((wt_cc*c_wt*((1+inf)^wt_l)) - wt_sv*c_wt))/ ( ((1+inf)^wt_l) - 1));

dcc_pv = ((inf*((pv_cc*c_pv*((1+inf)^pv_l)) - pv_sv*c_pv))/ ( ((1+inf)^pv_l) - 1));

doc_wt=((inf*wt_mc*c_wt)/ ( ((1+inf)^365) - 1));

doc_pv=((inf*pv_mc*c_pv)/ ( ((1+inf)^365) - 1));

z = z + dcc_wt + dcc_pv + doc_wt + doc_pv + dccc_bat ;
end

