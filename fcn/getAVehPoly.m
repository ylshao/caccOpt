function aVehPoly = getAVehPoly(vVeh, lambda2, lambda3, aoCoeff, boCoeff)
    
    V_OC = 201.6; %volt
    Q_BATT = 6.5*3600; % ampere*sec
    R_BATT = 0.003*6*28;  % ohm

    % get initial aVeh and pBatt
    % poly: s^3 - 6*s^2 - 11*s - 27 => [1 -6 -11 -27]
    pdAoPoly = [4*aoCoeff(12), 3*aoCoeff(11)*vVeh + 3*aoCoeff(9), ...
        2*aoCoeff(10)*vVeh^2 + 2*aoCoeff(8)*vVeh + 2*aoCoeff(6), ...
        aoCoeff(7)*vVeh^2 + aoCoeff(5)*vVeh + aoCoeff(3)];

    pdBoPoly = [2*boCoeff(6), boCoeff(5)*vVeh + boCoeff(3)]; 

    aoPoly = [aoCoeff(12), aoCoeff(11)*vVeh + aoCoeff(9), ...
        aoCoeff(10)*vVeh^2 + aoCoeff(8)*vVeh + aoCoeff(6), ...
        aoCoeff(7)*vVeh^2 + aoCoeff(5)*vVeh + aoCoeff(3), ...
        aoCoeff(4)*vVeh^2 + aoCoeff(2)*vVeh + aoCoeff(1)];

    boPoly = [boCoeff(6), boCoeff(5)*vVeh + boCoeff(3), ...
        boCoeff(4)*vVeh^2 + boCoeff(2)*vVeh + boCoeff(1)];

    % get the polynomial equation for aVeh calculation

    partOne = -V_OC^2/(4*R_BATT)*Q_BATT^2*conv(pdAoPoly, conv(aoPoly, aoPoly));
    partTwo = lambda3^2/(4*R_BATT)*pdAoPoly;
    partThree = -Q_BATT^2*conv(pdBoPoly, conv(aoPoly, aoPoly));
    partFour = -lambda2*Q_BATT^2*conv(aoPoly, aoPoly);

    ORDER = max(numel(partOne), max(numel(partTwo), ...
        max(numel(partThree), numel(partFour))));

    aVehPoly = [zeros(1, ORDER-numel(partOne)), partOne] + ...
            [zeros(1, ORDER-numel(partTwo)), partTwo] + ...
            [zeros(1, ORDER-numel(partThree)), partThree] + ...
            [zeros(1, ORDER-numel(partFour)), partFour];
%         
%          % get the initial acceleration
%     aVehPolyRoots = roots(aVehPoly);
%     aVehPolyRealRoots = aVehPolyRoots(aVehPolyRoots == real(aVehPolyRoots))';
%     
% %     % METHOD 1: choose the smallest one here
% %     [~, minAVehInd] = min(abs(aVehPolyRealRoots));
% %     aVeh(i+1) = aVehPolyRealRoots(minAVehInd); 
%     
% %    % METHOD 2: choose the one that closest to last aVeh
% %     [~, minAVehInd] = min(abs(aVehPolyRealRoots - aVeh(i)));
% %     aVeh(i+1) = aVehPolyRealRoots(minAVehInd); 
% 
% %     % METHOD 3: choose the one that closest to aPre
% %     [~, minAVehInd] = min(abs(aVehPolyRealRoots - aPre(i+1)));
% %     aVeh(i+1) = aVehPolyRealRoots(minAVehInd); 
% 
%     % METHOD 4: find correct aVeh through Hamiltonian
%     ROOTS_NUM = numel(aVehPolyRealRoots);
%     if ROOTS_NUM == 1
%         aVeh = aVehPolyRealRoots; 
%         pBatt = 1/(4*R_BATT)*(V_OC^2 - (LAMBDA3/(aoFit(vVeh, aVeh(i+1))*Q_BATT))^2);
%     else
%         hamilArray = nan(numel(ROOTS_NUM), 1);
%         for iRoots = 1:ROOTS_NUM
%             aVehCur = aVehPolyRealRoots(iRoots);
%             pBattCur = 1/(4*R_BATT)*(V_OC^2 - (LAMBDA3/(aoFit(vVeh, aVehCur)*Q_BATT))^2);
%             fuelCons = fuelConsFcn(vVeh, aVehCur, pBattCur);
%             stateOne = LAMBDA1*(vPre(i+1)-vVeh);
%             stateTwo = lambda2*aVehCur;
%             stateThree = LAMBDA3*(-V_OC+sqrt(V_OC^2-4*R_BATT*pBattCur))/(2*R_BATT*Q_BATT);
%             
%             hamilArray(iRoots) = fuelCons + stateOne + stateTwo + stateThree;
%         end
%         [~, minAVehInd] = min(hamilArray);
%         aVeh = aVehPolyRealRoots(minAVehInd); 
%     end
   
        
        
% %% Sanity Check
% load('fuelConsFit_150624')
% 
% aoVehOne = polyval(aoPoly, linspace(-3, 3, 20));
% aoVehTwo = aoFit(ones(1, 20)*vVeh, linspace(-3, 3, 20));
% if abs(min(aoVehOne - aoVehTwo)) > 1e-10
%     debug = 1;
% end
% 
% figure; hold on;
% plot(aoVehOne, 'b')
% plot(aoVehTwo, 'ro')
% 
% boVehOne = polyval(boPoly, linspace(-3, 3, 20));
% boVehTwo = boFit(ones(1, 20)*vVeh, linspace(-3, 3, 20));
% if abs(min(boVehOne - boVehTwo)) > 1e-10
%     debug = 1;
% end
% figure; hold on;
% plot(boVehOne, 'b')
% plot(boVehTwo, 'ro')
% 
% pdAoOne = polyval(pdAoPoly, linspace(-3, 3, 20));
% pdAoTwo = pdAoAVehFcn(ones(1, 20)*vVeh, linspace(-3, 3, 20));
% if abs(min(pdAoOne - pdAoTwo)) > 1e-10
%     debug = 1;
% end
% figure; hold on;
% plot(pdAoOne, 'b')
% plot(pdAoTwo, 'ro')
% 
% pdBoOne = polyval(pdBoPoly, linspace(-3, 3, 20));
% pdBoTwo = pdBoAVehFcn(ones(1, 20)*vVeh, linspace(-3, 3, 20));
% if abs(min(pdBoOne - pdBoTwo)) > 1e-10
%     debug = 1;
% end
% figure; hold on;
% plot(pdBoOne, 'b')
% plot(pdBoTwo, 'ro')