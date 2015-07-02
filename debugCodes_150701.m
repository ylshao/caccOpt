%% debug aoFit, boFit and fuelCons
vVehMet = linspace(0, 60, 301)*MPH_2_KMPH*1000/3600;
aVehMet = linspace(-15, 15, 601)*MPH_2_KMPH*1000/3600;

pBatt = 20160;
    [vVehPlot, aVehPlot] = meshgrid(vVehMet, aVehMet); % 

    thisAoPlot = aoFitFullFcn(reshape(vVehPlot, [], 1), reshape(aVehPlot, [], 1));
    thisBoPlot = boFitFullFcn(reshape(vVehPlot, [], 1), reshape(aVehPlot, [], 1));
   
%     fuelConsPlot = thisAoPlot*pBatt+thisBoPlot;
fuelConsPlot = fuelConsFullFcn(reshape(vVehPlot, [], 1), reshape(aVehPlot, [], 1), ones(301*601, 1)*pBatt);
%     fuelConsPlot(fuelConsPlot < 0) = 0;
    fuelConsPlot = reshape(fuelConsPlot, 601, 301);
figure; 
h = surf(vVehPlot, aVehPlot, fuelConsPlot);
set(h, 'LineStyle', 'none')
xlabel('vVeh [m/s]')
ylabel('aVeh [m/s^2]')
zlabel('fuel cons [g/s]')

a= 1