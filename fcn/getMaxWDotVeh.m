function wDotVehPlot = getMaxWDotVeh(tEngMap, wEngMap, Wv, Av)

TRAN_EFF = 0.9;
Mv = 1400; %kg
Rtire = 0.3107; %meter
Atire = 2.33; %m^2
Cd = 0.26;
rou= 1.202;
ftire = 0.00475;
g = 9.8; %N*m/sec
Jv = Mv*Rtire^2;
Voc = 201.6; %volt
Qbatt = 6.5*3600; % ampere*sec
Rbatt = 0.003*6*28;  % ohm
nm = 0.85;
ng = 0.85;
phi =0;
S = 30;
Rr = 78;
R = Rr/S;
K = 4.113;

PBATT_MAX_COEFF = 0.5;

    pEngMin = tEngMap(1)*wEngMap(1);
    pEngMax = tEngMap(end)*wEngMap(end);
    pBattMax = Voc^2/(4*Rbatt)*PBATT_MAX_COEFF;
    pVehMax = (pEngMax + pBattMax*nm);
    WvResamp = linspace(min(Wv), max(Wv), 50)'; % [rad/s]
    AvResamp = linspace(min(Av), max(Av), 1000)'; % [rad/s^2]
    wDotVehMaxSpan = nan(numel(WvResamp), 1);

    POS_IND = find(AvResamp > 0, 1, 'first');
    for iWVeh = 1:numel(WvResamp)
        thisWVeh = WvResamp(iWVeh);
        if thisWVeh < 0.02
            thisWVeh = 0.02;
        end

        WDOT_MAX_IND = nan;
        for iWDotVeh = POS_IND:numel(AvResamp)
            thisWDotVeh = AvResamp(iWDotVeh);
            thisTVeh = (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*thisWVeh^2 + Mv*thisWDotVeh*Rtire)*Rtire/TRAN_EFF;
            thisPBattMin = getPbatt(tEngMap(end), wEngMap(end), thisTVeh, thisWVeh);
            if thisPBattMin > pBattMax
                WDOT_MAX_IND = iWDotVeh-1;
                break;
            end
        end    

        if isnan(WDOT_MAX_IND);
            wDotVehMaxSpan(iWVeh) = AvResamp(end);
        else
            wDotVehMaxSpan(iWVeh) = AvResamp(WDOT_MAX_IND-1);
        end
    %     thisWDotVeh = (pVehMax - ftire*Mv*g*cos(phi)*Rtire*thisWVeh - ...
    %         0.5*rou*Cd*Atire*thisWVeh^3*Rtire^3)/(Mv*Rtire^2*thisWVeh);
    %     
    %     % verify the result
    %     thisTVeh = (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*thisWVeh^2 + Mv*thisWDotVeh*Rtire)*Rtire;
    %     thisPVehMax = thisTVeh*thisWVeh;
    %     if abs(thisPVehMax - pVehMax) > 1e-5
    %         debug = 1;
    %     end
    %     wDotVehMaxSpan(iWVeh) = min(thisWDotVeh, max(wDotVehSpan));

    end
    
    wDotVehPlot = nan(numel(Av), numel(Wv));
    for i = 1:length(Wv)
        wDotVehMaxCur = interp1(WvResamp, wDotVehMaxSpan, Wv(i), 'pchip');
        wDotVehPlot(:, i) = linspace(-15*1.61*1000/3600/Rtire, wDotVehMaxCur, numel(Av))';
    end
%%
    vVehSpanKmph = WvResamp*Rtire*3600/1000; % [kMph]
    aVehMaxSpan = wDotVehMaxSpan*Rtire;
    figure;
    plot(vVehSpanKmph, aVehMaxSpan)
    xlabel('vVeh [kMph]')
    ylabel('aVeh [m/s^2]')
