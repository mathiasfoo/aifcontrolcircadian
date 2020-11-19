function dC = SB2004_AIFControl_ODE(t,C)


global theta Rf eta theta1 theta2 gammaC


dC = zeros(9,1);

% mPER2CRY
dC(1) = theta(1)*C(3)^(theta(2))*C(7)^(theta(3)) - theta(4)*C(1);

% pcPER2CRY
dC(2) = theta(5)*C(1)^(theta(6)) - theta(7)*C(2);

% pnPER2CRY
dC(3) = theta(8)*C(2)^(theta(9)) - theta(10)*C(3);

% mBMAL1
% Without AIF Control
% dC(4) = theta(11)*C(3)^(theta(12)*0) - theta(13)*C(4);
% With AIF Control
dC(4) = theta1*theta(11)*C(8)^(theta(12)) - theta(13)*C(4);

% pcBMAL1
dC(5) = theta(14)*C(4)^(theta(15)) - theta(16)*C(5);

% pnBMAL1
dC(6) = theta(17)*C(5)^(theta(18)) - theta(19)*C(6);

% paBMAL1
dC(7) = theta(20)*C(6)^(theta(21)) - theta(22)*C(7);

%% AIF Controller

% Z1
dC(8) = real(Rf - eta*C(8)*C(9) - gammaC*C(8));
% Z2
dC(9) = real(theta2*C(4) - eta*C(8)*C(9) - gammaC*C(9));
