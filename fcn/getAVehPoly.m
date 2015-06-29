function [aVehPolyRealRootsCell, aVehRoots] = getAVehPoly(vVeh, lambda2, lambda3, FitPara)
    
    V_OC = 201.6; %volt
    Q_BATT = 6.5*3600; % ampere*sec
    R_BATT = 0.003*6*28;  % ohm
    
AVEH_MAX = 5;
AVEH_MIN = -5;
    PIECE_ONE_LIM = -0.268333333333333;
    PIECE_TWO_LIM = 0.626111111111111;

aVehRoots = [];
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

            pdBoPoly = [2*boCoeff(6), boCoeff(5)*vVeh + boCoeff(3)]; 
            % ao = p1*a + p0
            aoPoly = [aoCoeff(9)*vVeh.^3 + aoCoeff(7)*vVeh.^2 + ...
                aoCoeff(5)*vVeh + aoCoeff(3), aoCoeff(8)*vVeh.^4 + ...
                aoCoeff(6)*vVeh.^3 + aoCoeff(4)*vVeh.^2 + ...
                aoCoeff(2)*vVeh + aoCoeff(1)];

%             boPoly = [boCoeff(6), boCoeff(5)*vVeh + boCoeff(3), ...
%                 boCoeff(4)*vVeh^2 + boCoeff(2)*vVeh + boCoeff(1)];
        else
            pdAoPoly = [4*aoCoeff(12), 3*aoCoeff(11)*vVeh + 3*aoCoeff(9), ...
                2*aoCoeff(10)*vVeh^2 + 2*aoCoeff(8)*vVeh + 2*aoCoeff(6), ...
                aoCoeff(7)*vVeh^2 + aoCoeff(5)*vVeh + aoCoeff(3)];

            pdBoPoly = [2*boCoeff(6), boCoeff(5)*vVeh + boCoeff(3)]; 
            % ao = p4*a^4 + p3*a^3 + p2*a^2+ p1*a + p0
            aoPoly = [aoCoeff(12), aoCoeff(11)*vVeh + aoCoeff(9), ...
                aoCoeff(10)*vVeh^2 + aoCoeff(8)*vVeh + aoCoeff(6), ...
                aoCoeff(7)*vVeh^2 + aoCoeff(5)*vVeh + aoCoeff(3), ...
                aoCoeff(4)*vVeh^2 + aoCoeff(2)*vVeh + aoCoeff(1)];

%             boPoly = [boCoeff(6), boCoeff(5)*vVeh + boCoeff(3), ...
%                 boCoeff(4)*vVeh^2 + boCoeff(2)*vVeh + boCoeff(1)];
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
        aVehPolyRealRoots = aVehPolyRoots(aVehPolyRoots == real(aVehPolyRoots))';
        
%         switch i
%             case 1
%             aVehPolyRealRoots = ...
%                 aVehPolyRealRoots(AVEH_MIN <= aVehPolyRealRoots & ...
%                 aVehPolyRealRoots < PIECE_ONE_LIM);
%             case 2
%             aVehPolyRealRoots = ...
%                 aVehPolyRealRoots(aVehPolyRealRoots >= PIECE_ONE_LIM & ...
%                 aVehPolyRealRoots <= PIECE_TWO_LIM);    
%             case 3
%             aVehPolyRealRoots = ...
%                 aVehPolyRealRoots(aVehPolyRealRoots > PIECE_TWO_LIM & ...
%                 aVehPolyRealRoots <= AVEH_MAX);
%         end
        if ~isempty(aVehPolyRealRoots)
            aVehPolyRealRootsCell{i} = aVehPolyRealRoots;
            aVehRoots = [aVehRoots aVehPolyRealRoots];
        end
    end