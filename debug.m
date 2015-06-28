
v = vVehCur;
a = 0.35;
pBatt = -2000;
lambda2 = -2.0195;
lambda3 = -316.8177;

aV = 0.3531;
pBatt = -2493.0715;
lambda2 = -2.0152;
lambda3 = -320.3229;

fuelConsDebug= fuelConsFcn(vVehCur, a, pBatt);
stateOneDebug = lambda1*(vPre-vVehCur);
stateTwoDebug = lambda2*a;
stateThreeDebug = lambda3*(-V_OC+sqrt(V_OC^2-4*R_BATT*pBatt))/(2*R_BATT*Q_BATT);

fprintf('fuelCons %8.4f, state 1 %8.4f, state 2 %8.4f, state 3 %8.4f\n', ...
    fuelConsDebug, stateOneDebug, stateTwoDebug, stateThreeDebug)

H = fuelConsDebug + stateOneDebug + stateTwoDebug + stateThreeDebug;
