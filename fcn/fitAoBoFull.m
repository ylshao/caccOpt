function FitPara = fitAoBoFull(mapData)
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
% fuelCons = ao*Pbatt + bo, if Pbatt is zero, bo should at least be larger
% than zero. Though seems like it should equal to the fuel cons @ minimal
% engine output
boListFit = boList;
MIN_BO = 0;
boListFit(boListFit < MIN_BO) = MIN_BO;
% convert to mph
vVehMph = wVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s]->[mph]
aVehMph = aVehList*R_TIRE/MPH_2_KMPH*3600/1000; % [rad/s^2]->[mph/s]

selAccelCell = [{1:145}, {145:165}, {165:numel(aVehList)}];
polyTypeAoCell = [{'poly41'}, {'poly24'}, {'poly41'}];
polyTypeBoCell = [{'poly32'}, {'poly32'}, {'poly32'}];
FitPara(numel(selAccelCell)) = struct;
for i = 1:numel(selAccelCell)
    selAccel = selAccelCell{i}; % -3mph-3mph, %%%165:deSamp:end, 101:201, 145:163 121:181
    aVehComfMph = aVehMph(selAccel); 
    aoSel = aoList(selAccel, :);
    boSel = boListFit(selAccel, :);

    % The selected region of ao bo mapping can use 
    %   mesh(vVehMapMph, aVehMapMph, aoMap) to plot
    [vVehMapMph, aVehMapMph] = meshgrid(vVehMph, aVehComfMph); % 
    aoMap = aoSel; %
    boMap = boSel; %
%     figure; mesh(vVehMapMph, aVehMapMph, aoMap)
    % for fitting, need to 1) reshape as vectors, 2) convert to metric units
    vVehFit = reshape(vVehMapMph*MPH_2_KMPH*1000/3600, [], 1); % [m/s]
    aVehFit = reshape(aVehMapMph*MPH_2_KMPH*1000/3600, [], 1); % [m/s^2]
    aoFit = reshape(aoMap, [], 1);
    boFit = reshape(boMap, [], 1);

    %% Fit ao bo using polynominals
    aoFitFcn = fit([vVehFit, aVehFit], aoFit, polyTypeAoCell{i});
    boFitFcn = fit([vVehFit, aVehFit], boFit, polyTypeBoCell{i});
    aoCoeff = coeffvalues(aoFitFcn);
    boCoeff = coeffvalues(boFitFcn);
    
    %% get rid of numerical errors
    NUMEL_ERROR_LIM = 1e-18;%
    for iAoInd = 1:numel(aoCoeff)
        if abs(aoCoeff(iAoInd)) < NUMEL_ERROR_LIM
            aoCoeff(iAoInd) = 0;
        end
    end
    NUMEL_ERROR_LIM = 1e-10;
    for iBoInd = 1:numel(boCoeff)
        if abs(boCoeff(iBoInd)) < NUMEL_ERROR_LIM
            boCoeff(iBoInd) = 0;
        end
    end
    %% save to output struct
    FitPara(i).aoFitFcn = aoFitFcn;
    FitPara(i).boFitFcn = boFitFcn;
    FitPara(i).aoCoeff = aoCoeff;
    FitPara(i).boCoeff = boCoeff;
    
    % plot
% vVehPlot = reshape(vVehFit, size(aoMap, 1), size(aoMap, 2));
% aVehPlot = reshape(aVehFit, size(aoMap, 1), size(aoMap, 2));
% aoFitPlot = aoFitFcn(vVehFit, aVehFit);
% aoFitPlot = reshape(aoFitPlot, size(aoMap, 1), size(aoMap, 2));
% figure; hold on
% plot3(vVehFit, aVehFit, aoFitFcn(vVehFit, aVehFit), 'x');
% mesh(vVehPlot, aVehPlot, reshape(aoFit, size(aoMap, 1), size(aoMap, 2)))
end

end