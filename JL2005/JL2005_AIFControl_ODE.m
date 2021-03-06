function dC = JL2005_AIFControl_ODE(t,C)


global EL theta Rf eta theta1 theta2 gammaC

dC = zeros(7,1);

% LHY mRNA
% Without AIF Control
% dC(1) = real(theta(1)*C(7)*EL + theta(2)*C(6)^(theta(3)*0) - theta(4)*C(1));
% With AIF Control
dC(1) = real(1*theta(1)*C(7)*EL + theta1*theta(2)*C(8)^(theta(3)) - theta(4)*C(1));

% LHY protein cytosol
dC(2) = real(theta(5)*C(1)^theta(6) - theta(7)*C(2));

% LHY protein nucleus
dC(3) = real(theta(8)*C(2)^theta(9) - theta(10)*C(3));

% TOC1 mRNA
dC(4) = real(theta(11)*C(3)^theta(12) - theta(13)*C(4));

% TOC1 protein cytosol
dC(5) = real(theta(14)*C(4)^theta(15) - theta(16)*C(5));

% TOC1 protein nucleus
dC(6) = real(theta(17)*C(5)^theta(18) - theta(19)*C(6));

% Protein P
dC(7) = real(theta(20)*(1 - EL) - theta(21)*C(7) - theta(22)*EL*C(7));

%% AIF Controller

% Z1
dC(8) = real(Rf - eta*C(8)*C(9) - gammaC*C(8));
% Z2
dC(9) = real(theta2*C(1) - eta*C(8)*C(9) - gammaC*C(9));
