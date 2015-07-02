function [hamil, hamilTraj] = getHamil(lambda1, lambda2, lambda3, lambda4, vPre, vVeh, aVeh, pBatt, distFollow, distFollowCstr, fuelConsFullFcn)

V_OC = 201.6; %volt
Q_BATT = 6.5*3600; % ampere*sec
R_BATT = 0.003*6*28;  % ohm
%% get initial Hamiltonian and hamiltonian trajectories
fuelCons = fuelConsFullFcn(vVeh, aVeh, pBatt);
stateOne = lambda1*(vPre-vVeh);
stateTwo = lambda2*aVeh;
stateThree = lambda3*(-V_OC+sqrt(V_OC^2-4*R_BATT*pBatt))/(2*R_BATT*Q_BATT);
stateFour = lambda4*distFollowCstr(distFollow);

hamil = fuelCons + stateOne + stateTwo + stateThree + stateFour;
hamilTraj = [fuelCons, stateOne, stateTwo, stateThree, stateFour, hamil];
