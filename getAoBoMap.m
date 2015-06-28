clear all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gasoline engine related map/model %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% engine fuel consumption model
wEngMap=[1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 4000]*2*pi/60;  % (rad/s), speed range of the engine
lbft2Nm=1.356; %conversion from lbft to Nm
tEngMap=[6.3 12.5 18.8 25.1 31.3 37.6 43.9 50.1 56.4 62.7 68.9 75.2]*lbft2Nm;  % (N*m), torque range of the engine

% (g/s), fuel use map indexed vertically by enginemap_spd and horizontally by enginemap_trq
fuelConsMap = [
 0.1513  0.1984  0.2455  0.2925  0.3396  0.3867  0.4338  0.4808  0.5279  0.5279  0.5279  0.5279 
 0.1834  0.2423  0.3011  0.3599  0.4188  0.4776  0.5365  0.5953  0.6541  0.6689  0.6689  0.6689 
 0.2145  0.2851  0.3557  0.4263  0.4969  0.5675  0.6381  0.7087  0.7793  0.8146  0.8146  0.8146 
 0.2451  0.3274  0.4098  0.4922  0.5746  0.6570  0.7393  0.8217  0.9041  0.9659  0.9659  0.9659 
 0.2759  0.3700  0.4642  0.5583  0.6525  0.7466  0.8408  0.9349  1.0291  1.1232  1.1232  1.1232 
 0.3076  0.4135  0.5194  0.6253  0.7312  0.8371  0.9430  1.0490  1.1549  1.2608  1.2873  1.2873 
 0.3407  0.4584  0.5761  0.6937  0.8114  0.9291  1.0468  1.1645  1.2822  1.3998  1.4587  1.4587 
 0.3773  0.5068  0.6362  0.7657  0.8951  1.0246  1.1540  1.2835  1.4129  1.5424  1.6395  1.6395 
 0.4200  0.5612  0.7024  0.8436  0.9849  1.1261  1.2673  1.4085  1.5497  1.6910  1.8322  1.8322 
 0.4701  0.6231  0.7761  0.9290  1.0820  1.2350  1.3880  1.5410  1.6940  1.8470  1.9999  2.0382 
 0.5290  0.6938  0.8585  1.0233  1.1880  1.3528  1.5175  1.6823  1.8470  2.0118  2.1766  2.2589 
 0.6789  0.8672  1.0555  1.2438  1.4321  1.6204  1.8087  1.9970  2.1852  2.3735  2.5618  2.7501 ];
[tEngPlot, wEngPlot]=meshgrid(tEngMap, wEngMap);
pEngKwPlot = tEngPlot.*wEngPlot/1000;
enginemapFuelConsGpkwh = fuelConsMap./pEngKwPlot*3600;

% Draw Max Tq Line
tEngMaxMap = [tEngMap(1) 77.2920 82.0380 84.7500 86.7840 89.3604 91.1232 92.8860 94.6488 96.4116 98.1744 99.9372 101.9712];

wEngMaxMap = [1000 1010 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 4000];

% Max Tq Line Data Resampling
MaxSp = 1000:1:4000;
MaxTq = interp1(wEngMaxMap,tEngMaxMap,MaxSp);
% from RPM to rad/s
MaxSp = MaxSp*2*pi/60;

% % % %%
% % % contourList = 0.25:0.1:2.45;
% % figure; clf; hold on
% % [C, h] = contour(tEng, wEng,enginemapFuelCons, 20);
% % set(h, 'ShowText', 'on')
% % % plot(MaxSp, MaxTq, 'g*')
% % % % [C, h] = contour(engineOmega,engineTorq,enginemapPowerKw, 20);
% % % % set(h, 'ShowText', 'on', 'LineColor', 'r')
% % % FONT_SIZE = 12;
% % % xlabel('We [RPM]', 'FontSize', FONT_SIZE)
% % % ylabel('Te [Nm]', 'FontSize', FONT_SIZE)
% % % title('Fuel-Rate Map [gram/sec]', 'FontSize', FONT_SIZE)
% % % set(gca, 'FontSize', FONT_SIZE)
% % % box on
% % %%
% % figure; clf; hold on
% % contourList = 0.25:0.1:2.45;
% % [C, h] = contour(enginemapOmega*60/(2*pi), enginemapTorq, enginemapFuelCons, contourList);
% % set(h, 'ShowText', 'on')
% % 
% % figure; surf(tEng, wEng*60/(2*pi), enginemapFuelCons)
%%
% Toyota Prius hybrid parameters              
Jec = (0.178+0.835+0.0062+0.0127); %kg*m^2
Jgs = 0.023; %kg*m^2
Jmr = 0.023; %kg*m^2
Kratio = 4.113;
%Mv = 1361; %kg
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
Hl = 42*1e6; % 42MJ/kg
phi =0;
S = 30;
R = 78;
% R = Rr/S;
K = 4.113;
kcom = [1 1 -1 -1;1 -1 1 -1];

SOCdiff = 0.01;

%
% vehVel 0-96
% vehAccel -6.7-6.7
MPH_2_KMPH = 1.60934;

% %%
% tVehFcn = @(vehOmega, vehOmegaDot) (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*vehOmega^2 + Mv*vehOmegaDot*Rtire)*Rtire;
% 
% tVehList = nan(numel(wDotVehList), numel(wVehList));
% pVehList = nan(numel(wDotVehList), numel(wVehList));
% for iW = 1:numel(wVehList)
%     for iWDot = 1:numel(wDotVehList)
%         tVehList(iWDot, iW) = tVehFcn(wVehList(iW), wDotVehList(iWDot));
%         pVehList(iWDot, iW) = tVehList(iWDot, iW)*wVehList(iW);
%     end
% end
% pVehListKw = pVehList/1000;
% [wVehPlot, wDotVehPlot] = meshgrid(wVehList, wDotVehList);
% %%
% figure(2)
% surf(wVehPlot, wDotVehPlot, pVehListKw)
% xlabel('wVeh [rad/s]')
% ylabel('wDotVeh [rad/s^2]');
% zlabel('pVeh [kW]')

% %%
% vehPowerMax = max(max(vehPowerList));
% vehPowerMin = min(min(vehPowerList));

%% Try to find appropriate accel and vel sets
TRAN_EFF = 0.9;
pEngMin = tEngMap(1)*wEngMap(1);
pEngMax = tEngMap(end)*wEngMap(end);
pBattMax = Voc^2/(4*Rbatt);
pVehMax = (pEngMax + pBattMax*nm);
wVehSpan = linspace(0, 60, 50)'*MPH_2_KMPH*1000/3600/Rtire; % [rad/s]
wDotVehSpan = linspace(0, 15, 1000)'*MPH_2_KMPH*1000/3600/Rtire; % [rad/s^2]
wDotVehMaxSpan = nan(numel(wVehSpan), 1);


for iWVeh = 1:numel(wVehSpan)
    thisWVeh = wVehSpan(iWVeh);
    if thisWVeh < 0.02
        thisWVeh = 0.02;
    end
    
    WDOT_MAX_IND = nan;
    for iWDotVeh = 2:numel(wDotVehSpan)
        thisWDotVeh = wDotVehSpan(iWDotVeh);
        thisTVeh = (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*thisWVeh^2 + Mv*thisWDotVeh*Rtire)*Rtire/TRAN_EFF;
        thisPBattMin = getPbatt(tEngMap(end), wEngMap(end), thisTVeh, thisWVeh);
        if thisPBattMin > pBattMax
            WDOT_MAX_IND = iWDotVeh-1;
            break;
        end
    end    
    
    if isnan(WDOT_MAX_IND);
        wDotVehMaxSpan(iWVeh) = wDotVehSpan(end);
    else
        wDotVehMaxSpan(iWVeh) = wDotVehSpan(WDOT_MAX_IND-1);
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
vVehSpanKmph = wVehSpan*Rtire*3600/1000; % [kMph]
aVehMaxSpan = wDotVehMaxSpan*Rtire;
figure;
plot(vVehSpanKmph, aVehMaxSpan)
xlabel('vVeh [kMph]')
ylabel('aVeh [m/s^2]')

% save('wDotVehMax_V2_150625', 'wVehSpan', 'wDotVehMaxSpan')
%% Iteration to get ao and bo
WDOT_VEH_LIM= 15*MPH_2_KMPH*1000/3600/Rtire; %[rad/s^2]
W_LIST_LEN = 5;
WDOT_LIST_LEN = 10;
wVehList = linspace(0, 60, 5)'*MPH_2_KMPH*1000/3600/Rtire; % [rad/s]
wDotVehList = linspace(-15, 15, 10)'*MPH_2_KMPH*1000/3600/Rtire; % [rad/s^2]
pBattList = linspace(-30, 10, 30)'*10^3; % [Watts]
wEngList = linspace(1000, 4000, 40)'*2*pi/60; % [rad/s]
aoList = nan(numel(wDotVehList), numel(wVehList)); 
boList = nan(numel(wDotVehList), numel(wVehList));
fprintf('---------------------------------------------------------\n')

wVehPlot = nan(numel(wDotVehList), numel(wVehList));
wDotVehPlot = nan(numel(wDotVehList), numel(wVehList));
for iW = 1:numel(wVehList)   
    wDotVehMaxCur = interp1(wVehSpan, wDotVehMaxSpan, wVehList(iW));
    wDotVehList = linspace(-WDOT_VEH_LIM, wDotVehMaxCur, WDOT_LIST_LEN);
    for iWDot = 1:numel(wDotVehList)

        % get current wVeh and current wDotVeh
        wVeh = wVehList(iW);
        wDotVeh = wDotVehList(iWDot);
        
        wVehPlot(iWDot, iW) = wVeh;
        wDotVehPlot(iWDot, iW) = wDotVeh;
        
        tVeh = (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*wVeh^2 + Mv*wDotVeh*Rtire)*Rtire*1/TRAN_EFF;
%         pVeh = pVehList(iWDot, iW);
        wMot = K*wVeh;
    
        if iW == 2 && iWDot == 7
            a = 1;
        end
        
        pBattCurMax = getPbatt(tEngMap(1), wEngMap(1), tVeh, wVeh);
        pBattCurMin = getPbatt(tEngMap(end), wEngMap(end), tVeh, wVeh);
        
%         if pBattCurMax > Voc^2/(4*Rbatt) && pBattCurMin < Voc^2/(4*Rbatt)
%             pBattCurMax = Voc^2/(4*Rbatt);
%             fprintf('high pBattCurMax\n')
%         end
        
                % display current wVeh and current wDotVeh
        fprintf('wVeh %d, wDotVeh %d, pBattMin %8.4f, pBattMax %8.4f, tReq %8.4f\n',...
            iW, iWDot, pBattCurMin, pBattCurMax, tVeh)
        pBattList = linspace(pBattCurMin, pBattCurMax, 30)';
        
        fuelConsEachPower = nan(numel(pBattList), 1);
        for iPBatt = 1:numel(pBattList)
            pBatt = pBattList(iPBatt);
            
            fuelConsEachOp = nan(numel(wEngList), 1);
            
            tEngList = nan(numel(wEngList), 1);
            for iWEng = 1:numel(wEngList)
                wEng = wEngList(iWEng);

                wGen = (wEng*(R+S)-wMot*R)/S;
                %  always negative
                %  if wGen >= 0 % genOmega*genTorq<0, as generator kg=1
                kg = sign(wGen);
                
                % find the correct sign of km
                km = -1;
                tEng = (nm^km*wMot*tVeh/K-pBatt)*(S+R)/...
                        (nm^km*wMot*R+ng^kg*wGen*S);
                tMot = tVeh/K - tEng*R/(S+R);
                if tMot*wMot < 0 
                    km = 1;
                    tEng = (nm^km*wMot*tVeh/K-pBatt)*(S+R)/...
                        (nm^km*wMot*R+ng^kg*wGen*S);
                    tMot = tVeh/K - tEng*R/(S+R);
                end
                
                tGen = -tEng*(S/(S+R));


                
                if wVeh > 0
                    debug = 1;
                end
                
                %
%                 if iW == 4 && iWDot == 3
                pEng = tEng*wEng;
                pVeh = tVeh*wVeh;
                pMot = tMot*wMot;
                pGen = tGen*wGen;
                
                pProd = pEng + pGen + pMot;
                if abs(pProd - pVeh) > 1e-5
                    debug = 1;
                end
                
%                 fprintf('pVeh %8.4f,  pEng %8.4f, pMot %8.4f, pGen %8.4f\n',...
%                     pVeh, pEng, pMot, pGen)
%                 fprintf('motor %d, generator %d\n', km, kg)
%                 
%                 end
                
                % save tEng, wEng curve
                tEngList(iWEng) = tEng;
                fuelConsEachOp(iWEng) = interp2(tEngMap, wEngMap, ...
                        fuelConsMap, tEng, wEng);
                %fprintf('%8.4f\n', tEng)
            end
            
            tEngResamp = interp1(wEngList,tEngList,MaxSp);
            wEngResamp = MaxSp;
            
            isValid = (MaxTq - tEngResamp) > 0 | abs(MaxTq - tEngResamp) < 1e-5;
            if isempty(isValid(isValid ~= 0)) % if all tEng is larger than MaxTq
                [minTqDiff, isValid] = min(abs(MaxTq - tEngResamp));
                minTqDiff
                if (minTqDiff > 1)
                    debug = 1;
                end
            end
            fuelConsResamp = interp2(tEngMap, wEngMap, ...
                fuelConsMap, tEngResamp(isValid), wEngResamp(isValid));
            if isempty(fuelConsResamp)
                debug = 1;
            end
            fuelConsEachPower(iPBatt) = nanmin(fuelConsResamp);
                        
%             fuelConsEachPower(iPBatt) = nanmin(fuelConsEachOp);
        end
        validInd = ~isnan(fuelConsEachPower);
        lineCoeff = polyfit(pBattList(validInd), fuelConsEachPower(validInd), 1);
        
        aoList(iWDot, iW) = lineCoeff(1);
        boList(iWDot, iW) = lineCoeff(2);
    end
end

%%
figure
surf(wVehPlot*Rtire*3600/1000, wDotVehPlot*Rtire, aoList*10^3)
xlabel('vVeh [kMph]')
ylabel('aVeh [m/s^2]');
zlabel('ao')

% %%
% figure;
% plot(pBattList, fuelConsEachPower, '*-')
% %%
% % Predicted Speed
% Wvp = dlmread('vehspeedpred.txt')*1.61/3.6/Rtire;
% 
% % Actual Speed
% Wva = dlmread('vehspeedact.txt')*1.61/3.6/Rtire;
% 
% %%
% figure; hold on
% plot(Wvp);
% plot(Wva, 'r')