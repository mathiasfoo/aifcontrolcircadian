function dC = JD2016_ESSModel_ODE(t,C)


global theta EL


dC = zeros(12,1);

% LHY mRNA
dC(1) = real(theta(1)*C(4)^(theta(2))*C(6)^(theta(3)) - theta(4)*C(1) + theta(5)*C(9)*EL);

% LHY protein
dC(2) = real(theta(6)*C(1) - theta(7)*C(2) + theta(8)*EL*C(1));

% P97 mRNA
dC(3) = real(theta(9)*C(2)^(theta(10))*C(6)^(theta(11))*C(8)^(theta(12)) - theta(13)*C(3) + theta(14)*C(9)*EL);

% P97 protein
dC(4) = real(theta(15)*C(3) - theta(16)*C(4) + theta(17)*EL*C(4));

% P51 mRNA
dC(5) = real(theta(18)*C(2)^(theta(19))*C(6)^(theta(20)) - theta(21)*C(5));
 
% P51 protein
dC(6) = real(theta(22)*C(5) - theta(23)*C(6) + theta(24)*EL*C(6));

% EL mRNA
dC(7) = real(theta(25)*C(2)^(theta(26))*C(6)^(theta(27))*C(8)^(theta(28))*EL - theta(29)*C(7));

% EL protein
dC(8) = real(theta(30)*C(7) - theta(31)*C(8) + theta(32)*EL*C(8));

% Protein P
dC(9) = real(theta(33)*(1-EL) - theta(34)*C(9)*(1-EL) + theta(35)*EL*C(9));

% PIF mRNA
dC(10) = real(theta(36)*C(8)^(theta(37)) - theta(38)*C(10));

% PIF protein
dC(11) = real(theta(39)*C(10) - theta(40)*C(11) + theta(41)*EL*C(11));

% HYP protein
dC(12) = real(theta(42)*C(11)^theta(43));

