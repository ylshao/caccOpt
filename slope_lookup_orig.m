format long
global Voc Rbatt Qbatt dt fo PbattMsv

nm = 0.85;
ng = 0.85;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gasoline engine related map/model %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% engine fuel consumption model
enginemap_spd=[1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 4000]*2*pi/60;  % (rad/s), speed range of the engine
lbft2Nm=1.356; %conversion from lbft to Nm
enginemap_trq=[6.3 12.5 18.8 25.1 31.3 37.6 43.9 50.1 56.4 62.7 68.9 75.2]*lbft2Nm;  % (N*m), torque range of the engine

% (g/s), fuel use map indexed vertically by enginemap_spd and horizontally by enginemap_trq
enginemap = [
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
[T,w]=meshgrid(enginemap_trq, enginemap_spd);
enginemap_kW=T.*w/1000;
enginemap_gpkWh = enginemap./enginemap_kW*3600;

% T = [zeros(length(enginemap_trq),1) T];
% w = [zeros(1,length(enginemap_trq)); w];
% enginemap = [zeros(1,length(enginemap_trq)); enginemap];
% enginemap = [zeros(length(enginemap_spd)+1,1) enginemap];

% % polynomial extrapolate enginemap Torque to 108 Nm
% [enginemap,enginemap_trq,enginemap_spd]= gridfit(T,w,enginemap,[enginemap_trq 108],enginemap_spd);
% enginemap_trq = enginemap_trq(1,:);
% enginemap_spd = enginemap_spd(:,1);

% Draw Max Tq Line
% MaxTq_pt = [enginemap_trq(1) 78 90 100 102 103 104 105 106 108 108]; % Nm
MaxTq_pt = [enginemap_trq(1) 77.2920 82.0380 84.7500 86.7840 89.3604 91.1232 92.8860 94.6488 96.4116 98.1744 99.9372 101.9712];

% MaxSp_pt = [1000 1001 1500 2000 2250 2500 2750 3000 3250 3500 4000]; % RPM
MaxSp_pt = [1000 1010 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 4000];

% Max Tq Line Data Resampling
MaxSp = 1000:1:4000;
MaxTq = interp1(MaxSp_pt,MaxTq_pt,MaxSp);
% from RPM to rad/s
MaxSp = MaxSp*2*pi/60;

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
Rr = 78;
R = Rr/S;
K = 4.113;
kcom = [1 1 -1 -1;1 -1 1 -1];

% % Sample Wv from 0 - 60 mph 
Wv = linspace(0,60,10)*1.61/3.6/Rtire;
% Sample Av from -15 - 15 mph/s 
Av = linspace(-15,15,10)*1.61/3.6/Rtire;

% % Sample Wv from 0 - 60 mph 
% Wv = linspace(0,60,601)*1.61/3.6/Rtire;
% % Sample Av from -15 - 15 mph/s 
% Av = linspace(-15,15,301)*1.61/3.6/Rtire;

dt = 1;

% We size
WeM = linspace(1000*2*pi/60, 4000*2*pi/60, 30);

% Pbatt iteration size
N = 30;

fM = [];
TM = [];
WM = [];
fMnew = [];
TMnew = [];
WMnew = [];
% fo = zeros(N,length(Wv));
% Teo = zeros(N,length(Wv));
% Weo = zeros(N,length(Wv));

countin = 0;
countout = 0;
c = 0;
d = 0;
cno = [];
cover = [];
Pbat_lim = [];
PbattMsv = [];
fo = [];
Teo = [];
Weo = [];
Mfo = [];
MPbattMsv = [];

for m = 1:length(Av);
    
    m
    Accreq = Rtire*Av(m);
    
    for i = 1:length(Wv);
        
        Wreq = Wv(i);
        
        if Wreq <= 0.02
            Treq = 0;
        else
            
            % During decel braking
            if Accreq < 0
                % Torque request felt by the motor (Forces have to be multiplied with Rtire)
                Treq1 = ( (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*Wreq^2 + Mv*Accreq)*Rtire )/K;
                
                % Max Torque available by the motor (50 kWatts)
                Tmot1 = -5e4/(K*Wreq);
                
                % Cap motor torque to 400 Nm in the low-speed regions
                if Tmot1 < -400
                    Tmot1 = -400;
                else
                    dummy =0;
                end
                
                % If requested torque is bigger than the motor max Torque
                if Treq1 < Tmot1
                    % torque before the ratio
                    Treq = Tmot1*K;
                else
                    % torque before the ratio
                    Treq = Treq1*K;
                end
                
            else
                % Torque request (Forces have to be multiplied with Rtire)
                Treq = (ftire*Mv*g*cos(phi) + 0.5*rou*Cd*Atire*Rtire^2*Wreq^2 + Mv*Accreq)*Rtire;
            end
            
        end
        
        % Transmission Efficiency = 90%
        Treq = Treq*1.11111;
        
        % ========================
        % Initialize Pbatt lo & hi
        % ========================
        
        % take Teng lower than enginemap_trq(1) so that Pbatt can cover left
        % bottom corner of engine map
        Teng = enginemap_trq(1);
        Weng = enginemap_spd(1);
        T1 = -Teng/(1+R);
        T2 = Treq/K - R/(1+R)*Teng;
        W1 = -R*K*Wreq + (1+R)*Weng;
        W2 = K*Wreq;
        
        if T1*W1 <= 0 %generator
            k1 = 1;
        elseif T1*W1 > 0 %motor
            k1 = -1;
        end
        
        if T2*W2 <= 0 %generator
            k2 = 1;
        elseif T2*W2 > 0 %motor
            k2 = -1;
        end
        
        Pbathi = ng^k1*T1*W1 + nm^k2*T2*W2;
        
        Teng = enginemap_trq(end);
        Weng = enginemap_spd(end);
        T1 = -Teng/(1+R);
        T2 = Treq/K - R/(1+R)*Teng;
        W1 = -R*K*Wreq + (1+R)*Weng;
        W2 = K*Wreq;
        
        if T1*W1 <= 0 %generator
            k1 = 1;
        elseif T1*W1 > 0 %motor
            k1 = -1;
        end
        
        if T2*W2 <= 0 %generator
            k2 = 1;
        elseif T2*W2 > 0 %motor
            k2 = -1;
        end
        
        Pbatlo = ng^k1*T1*W1 + nm^k2*T2*W2;
        
        % PbattM = NOT INCLUDING ENGINE SHUTOFF FOR ITERATION PURPOSES
        PbattM = linspace(Pbatlo,Pbathi,N);
        % ============================
        % END Initialize Pbatt lo & hi
        % ============================
        
        
        % ======================================
        % Calculate Pbatt when Engine Shuts Off
        % ======================================
        Teng = 0;
        Weng = 0;
        T1 = -Teng/(1+R);
        T2 = Treq/K - R/(1+R)*Teng;
        W1 = -R*K*Wreq + (1+R)*Weng;
        W2 = K*Wreq;
        
        if T1*W1 <= 0
            k1 = 1;
        elseif T1*W1 > 0
            k1 = -1;
        end
        
        if T2*W2 <= 0
            k2 = 1;
        elseif T2*W2 > 0
            k2 = -1;
        end
        
        Pbat_shutoff = nm^k1*T1*W1 + ng^k2*T2*W2;
        % ==========================================
        % END Calculate Pbatt when Engine Shuts Off
        % ==========================================
        
        % Combine
        
        % Pbat_lim = [shutoff, lo]
        Pbat_lim(i,:) = [Pbat_shutoff Pbatlo];
        
        % PbattMsv = [lo..hi..shutoff]
        PbattMsv(:,i) = [PbattM Pbat_shutoff];
        
        %     abc = Voc^2 - 4*Rbatt*PbattM;
        %     dSOC(i,:) = -1/(2*Qbatt*Rbatt) * (Voc - sign(abc).*(abc).^(.5));
        
        for j = 1:length(PbattM)
            
            Pbatt = PbattM(j);
            
            for k = 1:length(WeM)
                
                We = WeM(k);
                Wm = K*Wreq;
                Wg = -R*K*Wreq + (1+R)*We;
                
                for l = 1:length(kcom)
                    
                    km = kcom(1,l);
                    kg = kcom(2,l);
                    mtx = [ng^kg*Wg nm^km*Wm 0; 0 1 R/(1+R); 1 0 1/(1+R)];
                    T = inv(mtx)*[Pbatt; Treq/K; 0];
                    
                    % verify this!
                    Tgo = T(1);
                    Tmo = T(2);
                    Teo = T(3);
                    
                    % find T*W that correlates the the right value of k (-1 or 1)
                    
                    % Problem 1 : more than 1 combination of km & kg that work
                    % Answer  1 : Check cover
                    
                    % Problem 2 : If no Tm & Tg fullfil any of the condition, it will remain the same as the Previous Tm & Tg
                    % Answer  2 : Check cno
                    
                    if (Tmo*Wm*km <= 0)&&(Tgo*Wg*kg <= 0)
                        Tm = Tmo;
                        Tg = Tgo;
                        Te = Teo;
                        countin = countin + 1;
                    else
                        countout = countout + 1;
                    end
                    
                    % Check if NO Combination works!
                    if countout == 4
                        c = c + 1;
                        cno(c,:) = [i j k];
                        
                        % Check if more than 1 Combinations work!
                    elseif countin > 1
                        d = d + 1;
                        cover(d,:) = [i j k];
                    end
                    
                end
                
                countin = 0;
                countout = 0;
                
                fM(i,j,k) = interp2(enginemap_trq,enginemap_spd,enginemap,Te,We)*dt;
                TM(i,j,k) = Te;
                WM(i,j,k) = We;
                
            end
            
            % Use Interp1 to limit the torque available for each Pbatt
            % according to the Max Torque line :
            
            % 1) resample Te based on
            WMm = squeeze(WM(i,j,:));
            TMm = squeeze(TM(i,j,:));
            Tresam = interp1(WMm,TMm,MaxSp);
            Wresam = MaxSp;
            
            % Find min difference bet Tq (intersection)
            dTints = abs(Tresam - MaxTq);
            [a b] = min(dTints);
            
            % Intersection point T, W & interpolate fuel
            Wints = Wresam(b);
            Tints = Tresam(b);
            fints = interp2(enginemap_trq,enginemap_spd,enginemap,Tints,Wints)*dt;
            
            % find We < Wints
            idxbel = find(WM(i,j,:) < Wints);
            % set fuel consumption for We < Wints to be 999
            fM(i,j,idxbel) = 999;
            
            % find We > Wints
            idxabv = find(WM(i,j,:) >= Wints);
            
            % Reshape selected data in 3D into 1D array
            WMbel = squeeze(WM(i,j,idxbel));
            WMabv = squeeze(WM(i,j,idxabv));
            TMbel = squeeze(TM(i,j,idxbel));
            TMabv = squeeze(TM(i,j,idxabv));
            fMbel = squeeze(fM(i,j,idxbel));
            fMabv = squeeze(fM(i,j,idxabv));
            
            % append data at maxTorque line in matrix
            WMnew(i,j,:) = [WMbel; Wints; WMabv];
            TMnew(i,j,:) = [TMbel; Tints; TMabv];
            fMnew(i,j,:) = [fMbel; fints; fMabv];
            
            % find min fuel consumption
            [fmin idx] = min(fMnew(i,j,:));
            
            fo(j,i)  = fMnew(i,j,idx);
            Teo(j,i) = TMnew(i,j,idx);
            Weo(j,i) = WMnew(i,j,idx);
            
        end
        
        % set first//last row of fo (optimal fuel consumption to be 0.1513//2.8785
        % because if Pbatt is slightly below the fuel map, then
        % mf for min Pbatt will be NaN (will interrupt optimization process)
        [rowf colf] = size(fo);
        fo(1,:) = 2.7501;
        fo(rowf,:) = 0.1513;
        
        % include zero for engine shutoff
        Mfo(m,i,:) = [fo(:,i); 0];
        MPbattMsv(m,i,:) = PbattMsv(:,i);
        
    end
    
end

% For Engine Shut-Off, fuel consumption == 0
fo = [fo; zeros(1,colf)];


%% Calculate Optimal Input PbattF
% ===============================

% Find Slopes of M_dot vs Pbatt at each time step
aoM = [];
boM = [];
for n = 1:length(Av)
    
    for p = 1:length(Wv)
        
        % linear curve fitting
        ao = polyfit(squeeze(MPbattMsv(n,p,:)), squeeze(Mfo(n,p,:)), 1);
        
        aoM(n,p) = ao(1);
        boM(n,p) = ao(2);
        
    end
    
end

aoMfig = aoM;
boMfig = boM;
WvMfig = Wv;
AvMfig = Av;

% Save Lookup Table into .mat file
save('Lookup.mat','aoMfig','boMfig','WvMfig','AvMfig','Rtire')

% Compare plots (plot one on top of the other)
figure(1); 
mesh(WvMfig/1.61*3.6*Rtire,AvMfig/1.61*3.6*Rtire,boMfig*1e3)

datax = WvMfig/1.61*3.6*Rtire;
datay = AvMfig/1.61*3.6*Rtire;
dataz = boMfig*1e3;
dataz1 = aoMfig*1e3;
 set(gca,'fontsize',14)
 
 set(gca,'XTick',[0 10 20 30 40 50 60],'YTick',[-15 -10 -5 0 5 10 15])
 
% xlabel('Vehicle Speed (mph)','fontsize',15)
% % h=get(gca,'xlabel');
% % set(h,'rotation',25)
% 
% ylabel('Vehicle Acceleration (mph/s)','fontsize',15)
% % h=get(gca,'ylabel');
% % set(h,'rotation',-18)
% 
zlabel('a_0-Slope ((g/s)/kWatt)','fontsize',15)
% % title('a-Slope Lookup Table','fontsize',15)
% axis([0 60 -15 15])

% xlabh = get(gca,'xlabel');
% set(xlabh,'Position',get(xlabh,'Position') + [0 0.5 0])
% 
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [0.5 0 0])


