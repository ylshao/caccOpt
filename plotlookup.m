load('E:\Azrin Files\032713 a-slope Lookup\Lookup2.mat')

figure 
% aoM is multiplied by 1e3 to convert to unit of (g/s)/kWatt from
% (g/s)/Watt
h = surf(WvMfig/1.61*3.6*Rtire,AvMfig/1.61*3.6*Rtire,aoMfig*1e3);
set(h,'LineStyle','none')
set(gca,'fontsize',14)
set(gca,'XTick',[0 10 20 30 40 50 60],'YTick',[-15 -10 -5 0 5 10 15])

xlabel('Vehicle Speed (mph)','fontsize',15)
ylabel('Vehicle Acceleration (mph/s)','fontsize',15)
zlabel('a_0-Slope ((g/s)/kWatt)','fontsize',15)
% axis([0 60 -15 15])
%%
x = WvMfig/1.61*3.6*Rtire;
y = AvMfig/1.61*3.6*Rtire;
z = aoMfig*1e3; 

%%
deSamp = 5;
xxx = x(1:deSamp:end);
yyy = y(140:165);
zzz = z(140:165, 1:deSamp:end);
%%
figure
[yy, xx] = meshgrid(x, y);   
h = surf(yy, xx, z);
set(h,'LineStyle','none')
% colormap jet

%%
wVehList = reshape(xx, [], 1);
aVehList = reshape(yy, [], 1);
a0VehList = reshape(z, [], 1);
%%
% meshgrid x y z
listNum = numel(wVehList);
zFit = [ones(listNum, 1), wVehList, aVehList, wVehList.^2, wVehList.*aVehList, ...
    aVehList.^2, wVehList.^3, wVehList.^2.*aVehList, ];
f = @(x, y) p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + ...
    p30*x^3 + p21*x^2*y+ p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + ...
    p22*x^2*y^2+ p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y + p32*x^3*y^2+ ...
    p23*x^2*y^3 + p14*x*y^4;


sf = fit([x, y], z, 'poly22');
% plot(sf, [x, y], z)
