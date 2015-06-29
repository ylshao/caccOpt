function pdAoAVeh = getPdAoAVeh(v, a, FitPara)

PIECE_ONE_LIM = -0.863641240210278;
PIECE_TWO_LIM = 2.015162893823980;

if a < PIECE_ONE_LIM
    aoCoeff = FitPara(1).aoCoeff;
    pdAoAVehFcn = @(v, a) aoCoeff(3) + aoCoeff(5)*v + aoCoeff(7)*v.^2 + ...
        aoCoeff(9)*v.^3;
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    aoCoeff = FitPara(2).aoCoeff;
    pdAoAVehFcn = @(v, a) aoCoeff(3) + aoCoeff(5)*v + 2*aoCoeff(6)*a + ...
        aoCoeff(7)*v.^2 + 2*aoCoeff(8)*v.*a + 3*aoCoeff(9)*a.^2 + ...
        2*aoCoeff(10)*v.^2.*a + 3*aoCoeff(11)*v.*a.^2 + 4*aoCoeff(12)*a.^3;
elseif a > PIECE_TWO_LIM
    aoCoeff = FitPara(3).aoCoeff;
    pdAoAVehFcn = @(v, a) aoCoeff(3) + aoCoeff(5)*v + aoCoeff(7)*v.^2 + ...
        aoCoeff(9)*v.^3;
end
pdAoAVeh = pdAoAVehFcn(v, a);