mapData = 'lookup_aobo_150623';
load(mapData)
MPH_2_KMPH = 1.60934;
R_TIRE = 0.3107;
%% Prepare data for fitting
% convert to mph
vVehMph = wVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s]->[mph]
aVehMph = aVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s^2]->[mph/s]
 
% because the ao-mapping is sort of piecewise, only select a region of
% acceleration for polynomial fitting. Otherwise need to have maybe 3
% different piecewise polynomials
selAccel = 1:145; % -3mph-3mph, %%%165:deSamp:end, 101:201, 145:163 121:181
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
aoFitFcn = fit([vVehFit, aVehFit], aoFit, 'poly41');
boFitFcn = fit([vVehFit, aVehFit], boFit, 'poly22');
aoCoeff = coeffvalues(aoFitFcn);
boCoeff = coeffvalues(boFitFcn);

%% plot fuel consumption map
pBatt = -20000;

selInd = 1:145; % acceleration -3m/s^2-3m/s^2
[vVehPlot, aVehPlot] = meshgrid(vVehMph*MPH_2_KMPH*1000/3600, aVehMph(selInd)*MPH_2_KMPH*1000/3600); % 

aoPlot = aoFitFcn(reshape(vVehPlot, [], 1), reshape(aVehPlot, [], 1));
boPlot = boFitFcn(reshape(vVehPlot, [], 1), reshape(aVehPlot, [], 1));
fuelConsPlot = aoPlot*pBatt+boPlot;
aoPlot = reshape(aoPlot, size(vVehPlot, 1), size(vVehPlot, 2));
boPlot = reshape(boPlot, size(vVehPlot, 1), size(vVehPlot, 2));
fuelConsPlot = reshape(fuelConsPlot, size(vVehPlot, 1), size(vVehPlot, 2));
figure; 
h = surf(vVehPlot, aVehPlot, fuelConsPlot);
set(h, 'LineStyle', 'none')
xlabel('vVeh [m/s]')
ylabel('aVeh [m/s^2]')
zlabel('fuel cons [g/s]')
% ylim([-3.5 3.5])