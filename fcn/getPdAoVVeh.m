function pdAoVVeh = getPdAoVVeh(v, a, FitPara)

PIECE_ONE_LIM = -0.863641240210278;
PIECE_TWO_LIM = 2.015162893823980;

if a < PIECE_ONE_LIM
    aoCoeff = FitPara(1).aoCoeff;
    pdAoVVehFcn = @(v, a) aoCoeff(2) + 2*aoCoeff(4)*v + aoCoeff(5)*a + ...
    3*aoCoeff(6)*v.^2 + 2*aoCoeff(7)*v.*a + 4*aoCoeff(8).v.^3 +...
    3*aoCoeff(9)*v.^2.*a;
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    aoCoeff = FitPara(2).aoCoeff;
    pdAoVVehFcn = @(v, a) aoCoeff(2) + 2*aoCoeff(4)*v + aoCoeff(5)*a + ...
    2*aoCoeff(7)*a.*v + aoCoeff(8)*a.^2 + 2*aoCoeff(10)*a.^2.*v + ...
    aoCoeff(11)*a.^3;
elseif a > PIECE_TWO_LIM
    aoCoeff = FitPara(3).aoCoeff;
    pdAoVVehFcn = @(v, a) aoCoeff(2) + 2*aoCoeff(4)*v + aoCoeff(5)*a + ...
    3*aoCoeff(6)*v.^2 + 2*aoCoeff(7)*v.*a + 4*aoCoeff(8).v.^3 +...
    3*aoCoeff(9)*v.^2.*a;
end
pdAoVVeh = pdAoVVehFcn(v, a);