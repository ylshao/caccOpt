function pdBoVVeh = getPdBoVVeh(v, a, FitPara)

    PIECE_ONE_LIM = -0.268333333333333;
    PIECE_TWO_LIM = 0.626111111111111;

if a < PIECE_ONE_LIM
    boCoeff = FitPara(1).boCoeff;
    pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a + ...
        3*boCoeff(7)*v.^2 + 2*boCoeff(8)*a.*v + boCoeff(9)*a.^2;
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    boCoeff = FitPara(2).boCoeff;
    pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a + ...
        3*boCoeff(7)*v.^2 + 2*boCoeff(8)*a.*v + boCoeff(9)*a.^2;
elseif a > PIECE_TWO_LIM
    boCoeff = FitPara(3).boCoeff;
    pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a + ...
        3*boCoeff(7)*v.^2 + 2*boCoeff(8)*a.*v + boCoeff(9)*a.^2;
end
pdBoVVeh = pdBoVVehFcn(v, a);