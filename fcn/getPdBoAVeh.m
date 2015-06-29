function pdBoAVeh = getPdBoAVeh(v, a, FitPara)

PIECE_ONE_LIM = -0.863641240210278;
PIECE_TWO_LIM = 2.015162893823980;

if a < PIECE_ONE_LIM
    boCoeff = FitPara(1).boCoeff;
    pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a; 
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    boCoeff = FitPara(2).boCoeff;
    pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a; 
elseif a > PIECE_TWO_LIM
    boCoeff = FitPara(3).boCoeff;
    pdBoAVehFcn = @(v, a) boCoeff(3) + boCoeff(5)*v + 2*boCoeff(6)*a; 
end
pdBoAVeh = pdBoAVehFcn(v, a);