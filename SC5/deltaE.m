function dc = deltaE(E,Charge,discount,CC)



Kw = 0.00015;
de = Kw * Charge;

NECL = 0.2*E/de;
id = discount / (365);

SV = 40.86*E;

dc = (id*( (CC*(1+id)^NECL) - SV))/ ( (1+id)^NECL - 1);


%% method 2 (unsued)
% 
% DOD_ini = 100-SOC_ini;DOD_final = 100-SOC_final;
% 
% alpha = -5440.35;
% beta = 1191.54;
% 
% H_ini = alpha * log(DOD_ini/100) + beta;
% H_final = alpha * log(DOD_final/100) + beta;
% 
% K_D_ini = (1./DOD_ini) .* (1-0.8.^(1./H_ini));
% K_D_final = (1./DOD_final) .* (1-0.8.^(1./H_final));
% 
% de = E .* abs(DOD_ini.*K_D_ini - DOD_final.*K_D_final);




% plot3(SOC_ini_,SOC_final_,DC,'.')
% xlabel('SOC_ini')
% ylabel('SOC_final')
% zlabel('DE cost')

end