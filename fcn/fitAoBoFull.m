function [aoFitFcn, boFitFcn, aoCoeff, boCoeff] = fitAoBoFull(mapData)
% [aoFitFcn, boFitFcn, aoCoeff, boCoeff] = fitAoBo(mapData)
%
% mapData is the .mat file that contains the mapping of ao, bo. The
% function loads this mapping data, and use polynomial to fit the data as
% analytical forms. 
%
% aoFitFcn, boFitFcn are fitobject, can be considered as function handles of the fitting, usage:
%   ao(v, a) = aoFitFcn(v, a)
%
% aoCoeff, boCoeff are the coefficients of the polynomials, in the form of
% 1-by-n vector.

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
selAccel = 121:181; % -3mph-3mph, %%%165:deSamp:end, 101:201, 145:163 121:181
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

end