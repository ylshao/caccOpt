function [aVehPolyValidRootsCell, aVehPolyRealRootsCell] = getAVehPoly(distHist, vVeh, aVehHist, lambda1, lambda2, lambda3, FitPara)
    
    V_OC = 201.6; %volt
    Q_BATT = 6.5*3600; % ampere*sec
    R_BATT = 0.003*6*28;  % ohm
    
% AVEH_MAX = 3;
% AVEH_MIN = -3;
    PIECE_ONE_LIM = -0.268333333333333;
    PIECE_TWO_LIM = 0.626111111111111;

%     aVehRoots = [];
    aVehPolyValidRootsCell = cell(3, 1);
    aVehPolyRealRootsCell = cell(3, 1);
    for i = 1:3
        aoCoeff = FitPara(i).aoCoeff;
        boCoeff = FitPara(i).boCoeff;
        
        % get initial aVeh and pBatt
        % poly: s^3 - 6*s^2 - 11*s - 27 => [1 -6 -11 -27]      
        if i == 1 || i == 3
            % pdAoAVeh = p0
            pdAoPoly = aoCoeff(9)*vVeh.^3 + aoCoeff(7)*vVeh.^2 + ...
                aoCoeff(5)*vVeh + aoCoeff(3);

            pdBoPoly = [2*boCoeff(9)*vVeh + 2*boCoeff(6), ...
                boCoeff(8)*vVeh.^2 + boCoeff(5)*vVeh + boCoeff(3)]; 
            % ao = p1*a + p0
            aoPoly = [aoCoeff(9)*vVeh.^3 + aoCoeff(7)*vVeh.^2 + ...
                aoCoeff(5)*vVeh + aoCoeff(3), aoCoeff(8)*vVeh.^4 + ...
                aoCoeff(6)*vVeh.^3 + aoCoeff(4)*vVeh.^2 + ...
                aoCoeff(2)*vVeh + aoCoeff(1)];
        else
            pdAoPoly = [4*aoCoeff(12), 3*aoCoeff(11)*vVeh + 3*aoCoeff(9), ...
                2*aoCoeff(10)*vVeh^2 + 2*aoCoeff(8)*vVeh + 2*aoCoeff(6), ...
                aoCoeff(7)*vVeh^2 + aoCoeff(5)*vVeh + aoCoeff(3)];

            pdBoPoly = [2*boCoeff(9)*vVeh + 2*boCoeff(6), ...
                boCoeff(8)*vVeh.^2 + boCoeff(5)*vVeh + boCoeff(3)]; 
            % ao = p4*a^4 + p3*a^3 + p2*a^2+ p1*a + p0
            aoPoly = [aoCoeff(12), aoCoeff(11)*vVeh + aoCoeff(9), ...
                aoCoeff(10)*vVeh^2 + aoCoeff(8)*vVeh + aoCoeff(6), ...
                aoCoeff(7)*vVeh^2 + aoCoeff(5)*vVeh + aoCoeff(3), ...
                aoCoeff(4)*vVeh^2 + aoCoeff(2)*vVeh + aoCoeff(1)];

        end
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

        %%
%         if 
%             debug = 1;
%         end
        aVehPolyRoots = roots(aVehPoly);
        aVehPolyRealRoots = aVehPolyRoots(aVehPolyRoots == real(aVehPolyRoots));
        
        switch i
            case 1
                if abs(partOne) < 1e-10 & abs(partTwo) < 1e-10 & abs(sum(aVehPoly)) < 1e-10
                    aVehPolyRealRoots = aVehHist(2) + sign(aVehHist(2)-aVehHist(1))*0.1;
                end
                aVehPolyValidRoots = ...
                    aVehPolyRealRoots(aVehPolyRealRoots < PIECE_ONE_LIM) & ...
                    aVehPolyRealRoots(aVehPolyRealRoots > 10);
            case 2
                aVehPolyValidRoots = ...
                    aVehPolyRealRoots(aVehPolyRealRoots >= PIECE_ONE_LIM & ...
                    aVehPolyRealRoots <= PIECE_TWO_LIM);    
            case 3
                if abs(partOne) < 1e-10 & abs(partTwo) < 1e-10 & abs(sum(aVehPoly)) < 1e-10
                    aoPoly = aoCoeff(8)*vVeh.^4 + aoCoeff(6)*vVeh.^3 +...
                        aoCoeff(4)*vVeh.^2 + aoCoeff(2)*vVeh + aoCoeff(1);
                    
                    pdAoVVehPoly = aoCoeff(2) + 2*aoCoeff(4)*vVeh + ...
                        3*aoCoeff(6)*vVeh.^2 + + 4*aoCoeff(8)*vVeh.^3;
                    pdBoVVehPoly = [2*boCoeff(8)*vVeh + boCoeff(5), ...
                        boCoeff(2) + 2*boCoeff(4)*vVeh + 3*boCoeff(7)*vVeh.^2];
                    
                    subpartOne = -V_OC^2/(4*R_BATT)*conv(pdAoVVehPoly, conv(aoPoly, aoPoly));
                    subpartTwo = lambda3^2/(4*R_BATT*Q_BATT^2)*pdAoVVehPoly;
                    subpartThree = -conv(pdBoVVehPoly, conv(aoPoly, aoPoly));
                    subpartFour = lambda1*conv(aoPoly, aoPoly);
                    subpartFive = conv([boCoeff(5) + 2*boCoeff(8)*vVeh, 0], conv(aoPoly, aoPoly));
                    
                    SUB_ORDER = max(numel(subpartOne), max(numel(subpartTwo), ...
                        max(numel(subpartThree), max(numel(subpartFour), numel(subpartFive)))));

                    aVehSubPoly = [zeros(1, SUB_ORDER-numel(subpartOne)), subpartOne] + ...
                            [zeros(1, SUB_ORDER-numel(subpartTwo)), subpartTwo] + ...
                            [zeros(1, SUB_ORDER-numel(subpartThree)), subpartThree] + ...
                            [zeros(1, SUB_ORDER-numel(subpartFour)), subpartFour] + ...
                            [zeros(1, SUB_ORDER-numel(subpartFive)), subpartFive];
                        
% %                         % syms lam3
% % lam3 = lambda3;
% % aoVal = aoCoeff(8)*vVeh.^4 + aoCoeff(6)*vVeh.^3 + aoCoeff(4)*vVeh.^2 + ...
% %     aoCoeff(2)*vVeh + aoCoeff(1);
% % coOne = (2*aoCoeff(4) + 6*aoCoeff(6)*vVeh + 12*aoCoeff(8)*vVeh.^2)*...
% %     (-V_OC^2/(4*R_BATT)+lam3^2/(4*R_BATT*Q_BATT^2*aoVal^2));
% % coTwo = (aoCoeff(2) + 2*aoCoeff(4)*vVeh + 3*aoCoeff(6)*vVeh.^2 + ...
% %     4*aoCoeff(8)*vVeh.^3)^2*(-lam3^2/(2*R_BATT*Q_BATT^2*aoVal^3));
% % coThree = -(2*boCoeff(4) + 6*boCoeff(7)*vVeh);
% % 
% % coSum = coOne + coTwo + coThree
% % 
% % % solve(coSum == 0, lam3)
                        
                        
                    if aVehHist(1) < PIECE_TWO_LIM  & aVehHist(2) < PIECE_TWO_LIM  
                        aVehPolyRealRoots = PIECE_TWO_LIM + 0.1 + sign(aVehHist(2)-aVehHist(1))*0.05;
                    else
                        aVehPolyRealRoots = aVehHist(2) + sign(aVehHist(2)-aVehHist(1))*0.1;
                    end
                    fprintf('Solve accel poly stucked\n')
                end
                aVehPolyValidRoots = ...
                    aVehPolyRealRoots(aVehPolyRealRoots > PIECE_TWO_LIM);
        end
        if ~isempty(aVehPolyRealRoots)
            aVehPolyRealRootsCell{i} = aVehPolyRealRoots;
        end
        if ~isempty(aVehPolyValidRoots)
            aVehPolyValidRootsCell{i} = aVehPolyValidRoots;
%             aVehRoots = [aVehRoots aVehPolyValidRoots];
        end
    end
    
aVehPolyValidRoots = cell2mat(aVehPolyValidRootsCell);
if numel(aVehPolyValidRoots) < 1
%     if aVehHist(1) < PIECE_TWO_LIM  & aVehHist(2) < PIECE_TWO_LIM  
%         aVehPolyValidRootsCell{3} = PIECE_TWO_LIM + 0.1 + sign(distHist(2)-distHist(1))*0.05;
%     else
        aVehPolyValidRootsCell{3} = aVehHist(2) + sign(distHist(2)-distHist(1))*0.1;
        if aVehPolyValidRootsCell{3} < PIECE_TWO_LIM
            aVehPolyValidRootsCell{3} = [];
        end
%         if aVehPolyValidRootsCell{}
%     end
end