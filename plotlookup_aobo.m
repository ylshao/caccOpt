
aVehList = nan(1, 301);
wVehList = nan(1, 601);
aoList = nan(301, 601);
boList = nan(301, 601);

curBegInd = 1;
for i = 0:8
    load(['Lookup', num2str(i)])
    thisListLen = numel(AvMfig);
    aVehList(curBegInd:curBegInd+thisListLen-1) = AvMfig;
    aoList(curBegInd:curBegInd+thisListLen-1, :) = aoMfig;
    boList(curBegInd:curBegInd+thisListLen-1, :) = boMfig;
    curBegInd = curBegInd+thisListLen;
end
wVehList = WvMfig;
%%
figure(5); clf;
FONT_SIZE = 27;
% aoM is multiplied by 1e3 to convert to unit of (g/s)/kWatt from
% (g/s)/Watt
h = surf(wVehList/1.61*3.6*Rtire,aVehList/1.61*3.6*Rtire,aoList*1e3);
set(h,'LineStyle','none')
set(gca,'fontsize',14)
set(gca,'XTick',[0 10 20 30 40 50 60],'YTick',[-15 -10 -5 0 5 10 15])
xlabel('Vehicle Speed (mph)','fontsize',FONT_SIZE)
ylabel('Vehicle Acceleration (mph/s)','fontsize',FONT_SIZE)
zlabel('a_0-Slope ((g/s)/kWatt)','fontsize',FONT_SIZE)
h = get(gca,'XLabel'); % Handle of the x label
set(h,'Rotation',13)
h = get(gca,'yLabel'); % Handle of the x label
set(h,'Rotation',-20)

%%
figure(6); clf
% aoM is multiplied by 1e3 to convert to unit of (g/s)/kWatt from
% (g/s)/Watt
h = surf(wVehList/1.61*3.6*Rtire,aVehList/1.61*3.6*Rtire,boList);
set(h,'LineStyle','none')
set(gca,'fontsize',14)
set(gca,'XTick',[0 10 20 30 40 50 60],'YTick',[-15 -10 -5 0 5 10 15])
xlabel('Vehicle Speed (mph)','fontsize',FONT_SIZE)
ylabel('Vehicle Acceleration (mph/s)','fontsize',FONT_SIZE)
zlabel('b_0-Slope ((g/s)/kWatt)','fontsize',FONT_SIZE)
h = get(gca,'XLabel'); % Handle of the x label
set(h,'Rotation',13)
h = get(gca,'yLabel'); % Handle of the x label
set(h,'Rotation',-20)