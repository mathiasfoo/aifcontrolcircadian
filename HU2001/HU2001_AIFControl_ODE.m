function dC = HU2001_AIFControl_ODE(t,C)


global theta Rf eta theta1 theta2 gammaC

dC = zeros(12,1);


% mPER
% Without AIF Control
% dC(1) = real(theta(1)*C(10)^(theta(2)*0)*C(6)^(theta(3)) - theta(4)*C(1));
% With AIF Control
dC(1) = real(theta1*theta(1)*C(11)^(theta(2))*C(6)^(theta(3)) - theta(4)*C(1));

% pPER
dC(2) = real(theta(5)*C(1)^(theta(6)) - theta(7)*C(2) - theta(8)*C(2)*C(4) + theta(9)*C(5));

% mTIM
dC(3) = real(theta(10)*C(10)^(theta(11))*C(6)^(theta(12)) - theta(13)*C(3));

% pTIM
dC(4) = real(theta(14)*C(3)^(theta(15)) - theta(16)*C(4) - theta(17)*C(2)*C(4) + theta(18)*C(5));

% pcPT
dC(5) = real(theta(19)*C(2)^(theta(20))*C(4)^(theta(21)) - theta(22)*C(5));

% pnPT
dC(6) = real(theta(23)*C(5)^(theta(24)) - theta(25)*C(6));

% mCLK
dC(7) = real(theta(26)*C(10)^(theta(27)*1)*C(6)^(abs(theta(28))) - theta(29)*C(7));

% pCLK
dC(8) = real(theta(30)*C(7)^(theta(31)) - theta(32)*C(8) + theta(33)*C(9));

% pcCC
dC(9) = real(theta(34)*C(8)^(theta(35)) - theta(36)*C(9));

% pnCC
dC(10) = real(theta(37)*C(9)^(theta(38)) - theta(39)*C(10));


%% AIF Controller

% Z1
dC(11) = real(Rf - eta*C(11)*C(12) - gammaC*C(11));
% Z2
dC(12) = real(theta2*C(1) - eta*C(11)*C(12) - gammaC*C(12));
