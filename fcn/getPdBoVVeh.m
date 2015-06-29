function pdBoVVeh = getPdBoVVeh(v, a, FitPara)

PIECE_ONE_LIM = -0.863641240210278;
PIECE_TWO_LIM = 2.015162893823980;

if a < PIECE_ONE_LIM
    boCoeff = FitPara(1).boCoeff;
    pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a;
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    boCoeff = FitPara(2).boCoeff;
    pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a;
elseif a > PIECE_TWO_LIM
    boCoeff = FitPara(3).boCoeff;
    pdBoVVehFcn = @(v, a) boCoeff(2) + 2*boCoeff(4)*v + boCoeff(5)*a;
end
pdBoVVeh = pdBoVVehFcn(v, a);