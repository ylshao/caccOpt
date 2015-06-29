function [aoFitFcn, boFitFcn, aoCoeff, boCoeff] = getCurAoBoMap(a, FitPara)

PIECE_ONE_LIM = -0.863641240210278;
PIECE_TWO_LIM = 2.015162893823980;

curPiece = nan;
if a < PIECE_ONE_LIM
    curPiece = 1;
elseif a >= PIECE_ONE_LIM && a <= PIECE_TWO_LIM
    curPiece = 2;
elseif a > PIECE_TWO_LIM
    curPiece = 3;
end

aoFitFcn = FitPara(curPiece).aoFitFcn;
boFitFcn = FitPara(curPiece).boFitFcn;
aoCoeff = FitPara(curPiece).aoCoeff;
boCoeff = FitPara(curPiece).boCoeff;