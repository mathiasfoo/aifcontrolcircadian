function dC = AD2015_ESSModel_ODE(t,C)


global theta

dC = zeros(9,1);

% mFRQ
dC(1) = real(theta(1)*C(6)^(theta(2)) - theta(3)*C(1));

% pcFRQ
dC(2) = real(theta(4)*C(1) - theta(5)*C(2));

% pnFRQ
dC(3) = real(theta(6)*C(2) - theta(7)*C(3) - theta(8)*C(3)*C(6) + theta(9)*C(7));

% mWC1
dC(4) = real(theta(10)*C(9)^(theta(11)) - theta(12)*C(4));

% pcWC1
dC(5) = real(theta(13)*C(4)^(theta(14))*C(2)^(theta(15)) - theta(16)*C(5));

% pnWC1
dC(6) = real(theta(17)*C(5) - theta(18)*C(6) - theta(19)*C(3)*C(6) + theta(20)*C(7));

% pFRQWC1
dC(7) = real(theta(21)*C(3)*C(6) - theta(22)*C(7));

% mCSP1
dC(8) = real(theta(23)*C(6)^(theta(24))*C(9)^(theta(25)) - theta(26)*C(8));

% % pCSP1
dC(9) = theta(27)*C(8) - theta(28)*C(9);
