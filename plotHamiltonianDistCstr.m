clear all
%%
load('lookup_aobo_150623')

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

DIST_MAX = 60;
DIST_MIN = 5;
%% Prepare data for fitting
% convert to mph
vVehMph = wVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s]->[mph]
aVehMph = aVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s^2]->[mph/s]
 
% because the ao-mapping is sort of piecewise, only select a region of
% acceleration for polynomial fitting. Otherwise need to have maybe 3
% different piecewise polynomials
selAccel = 145:165; % -3mph-3mph, %%%165:deSamp:end, 101:201, 145:163 121:181
aVehComfMph = aVehMph(selAccel); 
aoComf = aoList(selAccel, :);
boComf = boList(selAccel, :);

% The selected region of ao bo mapping can use 
%   mesh(vVehMap, aVehMap, aoMap) to plot
[vVehMapMph, aVehMapMph] = meshgrid(vVehMph, aVehComfMph); % 
aoMap = aoComf; %
boMap = boComf; %

% for fitting, need to 1) reshape as vectors, 2) convert to metric units
vVehFit = reshape(vVehMapMph*MPH_2_KMPH*1000/3600, [], 1); % [m/s]
aVehFit = reshape(aVehMapMph*MPH_2_KMPH*1000/3600, [], 1); % [m/s^2]
aoFit = reshape(aoMap, [], 1);
boFit = reshape(boMap, [], 1);

%% Fit ao bo using polynominals
aoFitFcn = fit([vVehFit, aVehFit], aoFit, 'poly24');
boFitFcn = fit([vVehFit, aVehFit], boFit, 'poly22');
aoCoeff = coeffvalues(aoFitFcn);
boCoeff = coeffvalues(boFitFcn);

%% Create fcn handles with vVeh, aVeh as arguments
% partial derivative of ao w.r.t. aVeh
pdAoAVehFcn = @(v, a) aoCoeff(3) + aoCoeff(5)*v + 2*aoCoeff(6)*a + ...
    aoCoeff(7)*v.^2 + 2*aoCoeff(8)*v.*a + 3*aoCoeff(9)*a.^2 + ...
    2*aoCoeff(10)*v.^2.*a + 3*aoCoeff(11)*v.*a.^2 + 4*aoCoeff(12)*a.^3;
% partial derivative of bo w.r.t. aVeh
pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a; 
% partial derivative of ao w.r.t. vVeh
pdAoVVehFcn = @(v, a) aoCoeff(2) + 2*aoCoeff(4)*v + aoCoeff(5)*a + ...
    2*aoCoeff(7)*a.*v + aoCoeff(8)*a.^2 + 2*aoCoeff(10)*a.^2.*v + ...
    aoCoeff(11) * a.^3; 
% partial derivative of bo w.r.t. vVeh
pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a;

tVehFcn = @(v, a)(F_TIRE*M_VEH*g*cos(phi) + ...
    0.5*ROU*C_D*A_TIRE*R_TIRE^2*v/R_TIRE^2 + M_VEH*a)*R_TIRE;

fuelConsFcn = @(v, a, pBatt) aoFitFcn(v, a)*pBatt + boFitFcn(v, a);

% heaviFcn = @(x, Xo) max((x-Xo), 0)/(x-Xo);

distFollowCstr = @(dist) (dist - DIST_MIN).^2.*heaviside(DIST_MIN - dist) + ...
    (DIST_MAX - dist).^2.*heaviside(dist - DIST_MAX);

pdDistFollowCstr = @(dist) 2*(dist - DIST_MIN).*heaviside(DIST_MIN - dist) - ...
    2*(DIST_MAX - dist).*heaviside(dist - DIST_MAX);
%% Estimate one set of lambda-2 and lambda-3
% because aVeh and pBatt have physical meanings, so currently I am
% estimating the optimal control to estimate the range of costates
vVehEst = 16.1018937388889; % [m/s]
distFollowEst = 50; % [m]
aVehEst = 0.3; %0.7454; % [m/s^2]
pBattEst = 5000; %-16414.0418; % [Watts]

% aVehEst = aVehOpt; %0.7454; % [m/s^2]
% pBattEst = pBattOpt; %-16414.0418; % [Watts]

% aVehEst = 0.3531;
% pBattEst = -2493.0715;

wVehEst = vVehEst/R_TIRE;
wDotVehEst = aVehEst/R_TIRE;
tVehEst = tVehFcn(vVehEst, aVehEst);


pBattMax = getPbatt(tEngMap(1), wEngMap(1), tVehEst, wVehEst);
pBattMin = getPbatt(tEngMap(end), wEngMap(end), tVehEst, wVehEst);

% saturate pBatt if it is not selected properly
if pBattEst < pBattMin
    pBattEst = pBattMin;
    fprintf('pBatt too small, range [%8.4f, %8.4f]\n', pBattMin, pBattMax)
elseif pBattEst > pBattMax
    pBattEst = pBattMax;
    fprintf('pBatt too large, range [%8.4f, %8.4f]\n', pBattMin, pBattMax)
end

aoEst = aoFitFcn(vVehEst, aVehEst);
boEst = boFitFcn(vVehEst, aVehEst);

lambda3 = sqrt(V_OC^2 - 4*R_BATT*pBattEst)*aoEst*Q_BATT;
lambda2 = -pdAoAVehFcn(vVehEst, aVehEst)*pBattEst - pdBoAVehFcn(vVehEst, aVehEst);
lambda1 = 0;
fuelConsEst = fuelConsFcn(vVehEst, aVehEst, pBattEst);
fprintf('lambda2 %8.4f, lambda3 %8.4f, fuelCons [g/s] %8.4f\n', ...
    lambda2, lambda3, fuelConsEst)
%% Estimate the shape of Hamiltonianwhos
vPre = vVehEst+2;
vVehCur = vVehEst;
distFollowCur = distFollowEst;
stateCur = [distFollowCur vVehCur 0.6]; % state [d, v, SOC]'

% iteration for Hamiltonian calculation
ASPAN_LEN = 100;
PSPAN_LEN = 50;
aSpan = linspace(min(aVehFit), max(aVehFit), ASPAN_LEN);

hamilPlot = nan(ASPAN_LEN, PSPAN_LEN);
aPlot = nan(ASPAN_LEN, PSPAN_LEN);
pBattPlot = nan(ASPAN_LEN, PSPAN_LEN);
    
aSpan = linspace(-3, 3, ASPAN_LEN);
pBattSpan = linspace(-40000, V_OC^2/(4*R_BATT), PSPAN_LEN);

for iAVeh = 1:numel(aSpan)
    aVehCur = aSpan(iAVeh);
    wVehCur = vVehCur/R_TIRE;
    wDotVehCur = aVehCur/R_TIRE;
    tVehEst = tVehFcn(vVehCur, aVehCur);
  
%     pBattMax = min(getPbatt(tEngMap(1), wEngMap(1), tVehEst, wVehEst)...
%         , V_OC^2/(4*R_BATT));
%     pBattMin = getPbatt(tEngMap(end), wEngMap(end), tVehEst, wVehEst);
%     
%     pBattSpan = linspace(pBattMin, pBattMax, PSPAN_LEN);
    for iPBatt = 1:numel(pBattSpan)
        pBattCur = pBattSpan(iPBatt);
        
        % calc Hamiltonian part by part to compare the value of each part
        fuelConsCur = fuelConsFcn(vVehCur, aVehCur, pBattCur);
        stateOne = lambda1*(vPre-vVehCur);
        stateTwo = lambda2*aVehCur;
        stateThree = lambda3*(-V_OC+sqrt(V_OC^2-4*R_BATT*pBattCur))/(2*R_BATT*Q_BATT);
        
        if (V_OC^2-4*R_BATT*pBattCur < 0)
            debug = 1;
        end
        
%         fprintf('fuelCons %8.4f, state 1 %8.4f, state 2 %8.4f, state 3 %8.4f\n', ...
%             fuelConsCur, stateOne, stateTwo, stateThree)
%         
        % save the current hamiltonian and control inputs
        hamilPlot(iAVeh, iPBatt) = fuelConsCur + stateOne + stateTwo + stateThree;
        aPlot(iAVeh, iPBatt) = aVehCur;
        pBattPlot(iAVeh, iPBatt) = pBattCur;
    end
end

%% Plot the Hamiltonian
% surface plot
hamilPlot = abs(hamilPlot);
figure(123); clf 
surf(aPlot, pBattPlot, hamilPlot)
hold on

% get the minimum hamiltonian for current states
[~, hamilMinInd] = min(hamilPlot);
[~, hamilMinCol] = min(min(hamilPlot));
hamilMinRow = hamilMinInd(hamilMinCol);
hamilMinVal = hamilPlot(hamilMinRow, hamilMinCol);

% get the optimal control for current states
aVehOpt = aPlot(hamilMinRow, hamilMinCol);
pBattOpt = pBattPlot(hamilMinRow, hamilMinCol);
scatter3(aVehOpt, pBattOpt, hamilMinVal, 'r')
xlabel('accel [m/s^2]')
ylabel('pBatt [Watts]')
zlabel('hamiltonian')

fuelConsOpt = fuelConsFcn(vVehCur, aVehOpt, pBattOpt);
fprintf('optimal control aVeh %8.4f, pBatt %8.4f\n', aVehOpt, pBattOpt)
fprintf('fuelCon from %8.4f to %8.4f\n', fuelConsEst, fuelConsOpt)
fprintf('optimal Hamiltonian %8.4f\n', hamilMinVal)
% % scatter plot
% figure; 
% aPlotVec = reshape(aPlot, [], 1);
% pBattPlotVec = reshape(pBattPlot, [], 1);
% hamilPlotVec = reshape(hamilPlot, [], 1);
% scatter3(aPlotVec, pBattPlotVec, hamilPlotVec)
% hold on
% [~, hamilMinInd] = min(hamilPlotVec);
% scatter3(aPlotVec(hamilMinInd), pBattPlotVec(hamilMinInd), hamilPlotVec(hamilMinInd), 'r')
% grid on
% 
% xlabel('accel [m/s^2]')
% ylabel('pBatt [Watts]')
% zlabel('hamiltonian')


%%
% save('fuelConsFit_150624', 'aoFit', 'boFit', 'aoCoeff', 'boCoeff', 'pdAoAVehFcn', 'pdBoAVehFcn', 'pdAoVVehFcn', 'pdBoVVehFcn','tVehFcn', 'fuelConsFcn')
