
clear all
% clc
%% Engine Map
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


%% Parameters and Unit Conversion Const
MPH_2_KMPH = 1.60934;
R_TIRE = 0.3107;
V_OC = 201.6; %volt
Q_BATT = 6.5*3600; % ampere*sec
R_BATT = 0.003*6*28;  % ohm
K = 4.113;
Mv = 1361; %kg
M_VEH = 1400; %kg
A_TIRE = 2.33; %m^2
C_D = 0.26;
ROU= 1.202;
F_TIRE = 0.00475;
g = 9.8; %N*m/sec
phi = 0;

%% Get preceding vehicle speed
vPreList = load('vehspeedact.txt');
vPreList = vPreList*MPH_2_KMPH*1000/3600; % [m/s]
TS_PRE = 0.1;
timeList = 0:TS_PRE:(numel(vPreList)-1)*TS_PRE;

%% get fuel consumption fitting
load('fuelConsFit_150624')

%% 

% initial conditions
vPreInit = vPreList(1);
aPreInit = (vPreList(2) - vPreList(1))/TS_PRE;
statesInit = [10, max(vPreInit-2, 0), 0.6]'; % [d, v, SOC]
distInit = statesInit(1);
vVehInit = statesInit(2);
socInit = statesInit(3);
LAMBDA1 = 0;
lambda2Init = -1.46297308598051; % bizarre lambda value to get smooth initial acceleration 
LAMBDA3 = -227.036499568364; 

METHOD = 4;
    
% get polynomial for solving aVeh
aVehPoly = getAVehPoly(vVehInit, lambda2Init, LAMBDA3, aoCoeff, boCoeff);
    
%% get the initial acceleration
aVehPolyRoots = roots(aVehPoly);
aVehPolyRealRoots = aVehPolyRoots(aVehPolyRoots == real(aVehPolyRoots))';
[~, minAVehInd] = min(abs(aVehPolyRealRoots - aPreInit));
aVehInit = aVehPolyRealRoots(minAVehInd); % choose the smallest one here
pBattInit = 1/(4*R_BATT)*(V_OC^2 - (LAMBDA3/(aoFit(vVehInit, aVehInit)*Q_BATT))^2);
fprintf('aVehInit %8.4f, pBattInit %8.4f\n', aVehInit, pBattInit)

wVehInit = vVehInit/R_TIRE;
wDotVehInit = aVehInit/R_TIRE;
tVehInit = tVehFcn(vVehInit, aVehInit);
%% start simulation of the preceding vehicle performance
tic;
DT = 0.5;
SIMU_STEPS = floor(timeList(end)/DT);
STOP_FLAG = 0;

% constraints
AVEH_MAX = 3;
AVEH_MIN = -3;

PBATT_MAX = V_OC^2/(4*R_BATT);
PBATT_MIN = -40000;

% preceding vehicle
timeSimu = (0:DT:(SIMU_STEPS)*DT)';
vPre = [vPreInit; interp1(timeList, vPreList, timeSimu(2:end))];
distPre = cumsum(vPre*DT);
aPre = [nan; diff(vPre)/DT];

% initialize states store arrays
distFollow = [distInit; nan(SIMU_STEPS, 1)];
vVeh = [vVehInit; nan(SIMU_STEPS, 1)];
soc = [socInit; nan(SIMU_STEPS, 1)];
lambda2 = [lambda2Init; nan(SIMU_STEPS, 1)];
aVeh = [aVehInit; nan(SIMU_STEPS, 1)];
pBatt = [pBattInit; nan(SIMU_STEPS, 1)];
wVeh = [wVehInit; nan(SIMU_STEPS, 1)];
wDotVeh = [wDotVehInit; nan(SIMU_STEPS, 1)];
tVeh = [tVehInit; nan(SIMU_STEPS, 1)];

for i = 1:SIMU_STEPS
    vPre(i) = vPre(i);
    distFollow(i+1) = distFollow(i) + DT*(vPre(i) - vVeh(i));
    vVeh(i+1) = vVeh(i) + DT*aVeh(i);
    soc(i+1) = soc(i) + DT*(-V_OC+sqrt(V_OC^2-4*R_BATT*pBatt(i)))/(2*R_BATT*Q_BATT);
    lambda2(i+1) = lambda2(i) + DT*(-pdAoVVehFcn(vVeh(i), aVeh(i))*pBatt(i) - ...
                pdBoVVehFcn(vVeh(i), aVeh(i)) + LAMBDA1);
            
    aVehPoly = getAVehPoly(vVeh(i+1), lambda2(i+1), LAMBDA3, aoCoeff, boCoeff);

    % get the initial acceleration
    aVehPolyRoots = roots(aVehPoly);
    aVehPolyRealRoots = aVehPolyRoots(aVehPolyRoots == real(aVehPolyRoots))';
    
    switch METHOD 
        case 1
            % METHOD 1: choose the smallest one here
            [~, minAVehInd] = min(abs(aVehPolyRealRoots));
            aVeh(i+1) = aVehPolyRealRoots(minAVehInd);             
        case 2
            % METHOD 2: choose the one that closest to last aVeh
            [~, minAVehInd] = min(abs(aVehPolyRealRoots - aVeh(i)));
            aVeh(i+1) = aVehPolyRealRoots(minAVehInd);             
        case 3
            % METHOD 3: choose the one that closest to aPre
            [~, minAVehInd] = min(abs(aVehPolyRealRoots - aPre(i+1)));
            aVeh(i+1) = aVehPolyRealRoots(minAVehInd);             
        case 4
            % METHOD 4: find correct aVeh through Hamiltonian
            ROOTS_NUM = numel(aVehPolyRealRoots);
            if ROOTS_NUM == 1
                aVeh(i+1) = aVehPolyRealRoots; 
            else
                hamilArray = nan(numel(ROOTS_NUM), 1);
                for iRoots = 1:ROOTS_NUM
                    aVehCur = aVehPolyRealRoots(iRoots);
                    pBattCur = 1/(4*R_BATT)*(V_OC^2 - (LAMBDA3/(aoFit(vVeh(i+1), aVehCur)*Q_BATT))^2);
                    fuelCons = fuelConsFcn(vVeh(i+1), aVehCur, pBattCur);
                    stateOne = LAMBDA1*(vPre(i+1)-vVeh(i+1));
                    stateTwo = lambda2(i+1)*aVehCur;
                    stateThree = LAMBDA3*(-V_OC+sqrt(V_OC^2-4*R_BATT*pBattCur))/(2*R_BATT*Q_BATT);

                    hamilArray(iRoots) = fuelCons + stateOne + stateTwo + stateThree;
                end
                [~, minAVehInd] = min(hamilArray);
                aVeh(i+1) = aVehPolyRealRoots(minAVehInd); 
            end
            
    end

    pBatt(i+1) = 1/(4*R_BATT)*(V_OC^2 - (LAMBDA3/(aoFit(vVeh(i+1), aVeh(i+1))*Q_BATT))^2);

    wVeh(i+1) = vVeh(i+1)/R_TIRE;
    wDotVeh(i+1) = aVeh(i+1)/R_TIRE;
    tVeh(i+1) = tVehFcn(vVeh(i+1), aVeh(i+1));
    pBattCurMax = max(min(getPbatt(tEngMap(1), wEngMap(1), tVeh(i+1), wVeh(i+1)), PBATT_MAX), PBATT_MIN);
    pBattCurMin = min(max(getPbatt(tEngMap(end), wEngMap(end), tVeh(i+1), wVeh(i+1)), PBATT_MIN), PBATT_MAX);
    % simple saturation
    aVeh(i+1) = min(max(aVeh(i+1), AVEH_MIN), AVEH_MAX);
    pBatt(i+1) = min(max(pBatt(i+1), pBattCurMin), pBattCurMax);
    
    if pBatt(i+1) < -40000
        debug = 1;
    end
end
toc;

%% plotting

SUB_ALIGN = '23';
ax = [];
subNum = 1;
figure; 
ax(subNum) = subplot([SUB_ALIGN, num2str(subNum)]); subNum = subNum + 1;
plot(timeSimu, distFollow)
grid on
xlabel('time [sec]')
ylabel('following dist [m]')

ax(subNum) = subplot([SUB_ALIGN, num2str(subNum)]); subNum = subNum + 1;
plot(timeSimu, [vPre, vVeh])
legend('vPre', 'vVeh')
grid on
xlabel('time [sec]')
ylabel('velocity [m/s]')

ax(subNum) = subplot([SUB_ALIGN, num2str(subNum)]); subNum = subNum + 1;
plot(timeSimu, [aPre, aVeh])
legend('aPre', 'aVeh')
grid on
xlabel('time [sec]')
ylabel('acceleration [m/s^2]')

ax(subNum) = subplot([SUB_ALIGN, num2str(subNum)]); subNum = subNum + 1;
plot(timeSimu, pBatt/10^3)
grid on
xlabel('time [sec]')
ylabel('pBatt [kWatts]')

ax(subNum) = subplot([SUB_ALIGN, num2str(subNum)]); subNum = subNum + 1;
plot(timeSimu, soc)
grid on
xlabel('time [sec]')
ylabel('soc')

ax(subNum) = subplot([SUB_ALIGN, num2str(subNum)]); subNum = subNum + 1;
plot(timeSimu, lambda2)
grid on
xlabel('time [sec]')
ylabel('lambda2')

linkaxes(ax, 'x')