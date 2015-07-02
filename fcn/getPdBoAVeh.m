function pdBoAVeh = getPdBoAVeh(v, a, FitPara)

    PIECE_ONE_LIM = -0.268333333333333;
    PIECE_TWO_LIM = 0.626111111111111;

if a < PIECE_ONE_LIM
    boCoeff = FitPara(1).boCoeff;
    pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a + ...
        boCoeff(8)*v.^2 + 2*boCoeff(9)*v.*a; 
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    boCoeff = FitPara(2).boCoeff;
    pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a + ...
        boCoeff(8)*v.^2 + 2*boCoeff(9)*v.*a; 
elseif a > PIECE_TWO_LIM
    boCoeff = FitPara(3).boCoeff;
    pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a + ...
        boCoeff(8)*v.^2 + 2*boCoeff(9)*v.*a; 
end
pdBoAVeh = pdBoAVehFcn(v, a);