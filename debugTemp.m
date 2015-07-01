mapData = 'lookup_aobo_150623';
load(mapData)
MPH_2_KMPH = 1.60934;
R_TIRE = 0.3107;
%% Prepare data for fitting
% fuelCons = ao*Pbatt + bo, if Pbatt is zero, bo should at least be the
% lowest engine output
boListFit = boList;
boListFit(boListFit < 0.1513) = 0.1513;


% convert to mph
vVehMph = wVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s]->[mph]
aVehMph = aVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s^2]->[mph/s]
 
% because the ao-mapping is sort of piecewise, only select a region of
% acceleration for polynomial fitting. Otherwise need to have maybe 3
% different piecewise polynomials
selAccel = 145:165; % -3mph-3mph, %%%165:deSamp:end, 101:201, 145:163 121:181
aVehComfMph = aVehMph(selAccel); 
aoComf = aoList(selAccel, :);
boComf = boListFit(selAccel, :);

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
aoFitFcn = fit([vVehFit, aVehFit], aoFit, 'poly41');
boFitFcn = fit([vVehFit, aVehFit], boFit, 'poly22');
aoCoeff = coeffvalues(aoFitFcn);
boCoeff = coeffvalues(boFitFcn);

%% plot fuel consumption map
FitPara = fitAoBoFull(mapData);
pBatt = 20160; % -40000 - 20160

selInd = 1:301; % acceleration -3m/s^2-3m/s^2

selIndArray = [{1:145}, {145:165}, {165:301}];
% selIndArray = [{145:165}];
aoPlot = [];
boPlot = [];
fuelConsPlot = [];
aVehPlot = [];
vVehPlot = [];
vVehMet = vVehMph*MPH_2_KMPH*1000/3600;
aVehMet = aVehMph*MPH_2_KMPH*1000/3600;
for ind = 1:numel(selIndArray)
    thisInd = selIndArray{ind};
    [thisVVeh, thisAVeh] = meshgrid(vVehMet, aVehMet(thisInd)); % 
    aoFitFcn = FitPara(ind).aoFitFcn;
    boFitFcn = FitPara(ind).boFitFcn;
    thisAoPlot = aoFitFcn(reshape(thisVVeh, [], 1), reshape(thisAVeh, [], 1));
    thisBoPlot = boFitFcn(reshape(thisVVeh, [], 1), reshape(thisAVeh, [], 1));
   
%     thisAoPlot = [];
%     thisBoPlot = [];
%     for iVVeh = 1:numel(vVehMet)
%         for iAVeh = 1:numel(thisInd)
%             thisAoPlot = [thisAoPlot; aoList(thisInd(iAVeh), iVVeh)];
%             thisBoPlot = [thisBoPlot; boList(thisInd(iAVeh), iVVeh)];
%         end
%     end
    PLOT_ROW_NUM = numel(thisAoPlot)/601;
    PLOT_COL_NUM = 601;
    thisAoPlot = reshape(thisAoPlot, PLOT_ROW_NUM, PLOT_COL_NUM);
    thisBoPlot = reshape(thisBoPlot, PLOT_ROW_NUM, PLOT_COL_NUM);
    thisFuelConsPlot = thisAoPlot*pBatt+thisBoPlot;
    thisFuelConsPlot = reshape(thisFuelConsPlot, PLOT_ROW_NUM, PLOT_COL_NUM);
    
    aVehPlot = [aVehPlot; thisAVeh];
    vVehPlot = [vVehPlot; thisVVeh];
    aoPlot = [aoPlot; thisAoPlot];
    boPlot = [boPlot; thisBoPlot];
    fuelConsPlot = [fuelConsPlot; thisFuelConsPlot];
end
fuelConsPlot(fuelConsPlot < 0) = 0;
figure; 
h = surf(vVehPlot, aVehPlot, fuelConsPlot);
set(h, 'LineStyle', 'none')
xlabel('vVeh [m/s]')
ylabel('aVeh [m/s^2]')
zlabel('fuel cons [g/s]')
% ylim([-3.5 3.5])

%%
aoArray = reshape(aoPlot, [], 1);
boArray = reshape(boPlot, [], 1);
fuelConsArray = reshape(fuelConsPlot, [], 1);
aoAdjArray = (fuelConsArray - boArray)/pBatt;
figure;
plot([aoArray aoAdjArray], 'x')
legend('ao', 'aoAdj')
aoAdjPlot = reshape(aoAdjArray, 301, 601);

figure;
h = surf(vVehPlot, aVehPlot, aoAdjPlot*pBatt+boPlot);
set(h, 'LineStyle', 'none')
xlabel('vVeh [m/s]')
ylabel('aVeh [m/s^2]')
zlabel('fuel cons [g/s]')
%%
figure; 
h = mesh(vVehPlot, aVehPlot, boPlot);
% set(h, 'LineStyle', 'none')