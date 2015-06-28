function pBatt = getPbatt(tEng, wEng, tVeh, wVeh)
% 
% Teng = enginemap_trq(1);
%         Weng = enginemap_spd(1);
S = 30;
Rr = 78;
R = Rr/S;
K = 4.113;
nm = 0.85;
ng = 0.85;

        tGen = -tEng/(1+R);
        tMot = tVeh/K - R/(1+R)*tEng;
        wGen = -R*K*wVeh + (1+R)*wEng;
        wMot = K*wVeh;
        
        if tGen*wGen <= 0 %generator
            kg = 1;
        elseif tGen*wGen > 0 %motor
            kg = -1;
        end
        
        if tMot*wMot <= 0 %generator
            km = 1;
        elseif tMot*wMot > 0 %motor
            km = -1;
        end
        
        pBatt = ng^kg*tGen*wGen + nm^km*tMot*wMot;