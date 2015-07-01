function fuelCons = getFuelCons(v, a, pBatt, FitPara)

PIECE_ONE_LIM = -0.268333333333333;
PIECE_TWO_LIM = 0.626111111111111;

if a < PIECE_ONE_LIM
    aoFitFcn = FitPara(1).aoFitFcn;
    boFitFcn = FitPara(1).boFitFcn;
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    aoFitFcn = FitPara(2).aoFitFcn;
    boFitFcn = FitPara(2).boFitFcn;
elseif a > PIECE_TWO_LIM
    aoFitFcn = FitPara(3).aoFitFcn;
    boFitFcn = FitPara(3).boFitFcn;
end
% fuel consumption handle
fuelCons = max(aoFitFcn(v, a)*pBatt + boFitFcn(v, a), 0);