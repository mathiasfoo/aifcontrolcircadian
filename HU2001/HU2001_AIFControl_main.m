clc
clear all
close all

%% Loading all the experimental WT data under different light condition.


global theta Rf eta theta1 theta2 gammaC

load drosophila2001LLdata.mat
load drosophila2001referencev1.mat 

eta = 105.6; 
theta1 = 5;
theta2 = 1;
gammaC = 1.5;

%% Model Parameters
prmPER = [0.4618,1.2751,-0.013493,0.5238];
prpPER = [0.43892,1.003,0.70352,1.3918,1.7507];
prmTIM = [0.37656,1.4944,-0.019446,0.48556];
prpTIM = [0.34403,1.0101,0.52987,1.3262,1.6196];
prpcPT = [1.441,1.4995,0.86091,2.1687];
prpnPT = [0.2821,1.2101,0.24144];
prmCLK = [0.29651,-1.7563,0.0023817,0.15121];
prpCLK = [0.17683,1.6382,2.2927,2.3236];
prpcCC = [1.4296,1.2849,2.5544];
prpnCC = [0.17936,1.7833,0.16931];

nonlineartheta = [prmPER prpPER prmTIM prpTIM prpcPT prpnPT prmCLK prpCLK prpcCC prpnCC];

theta = nonlineartheta;
GeneProteinLevelFull = [];

% return

%% Initial condition

C = [mPER(1) pPER(1) mTIM(1) pTIM(1) pcPT(1) pnPT(1) mCLK(1) pCLK(1) pcCC(1) pnCC(1) 0*ones(1,2)];
Cinit = [mPER(1) pPER(1) mTIM(1) pTIM(1) pcPT(1) pnPT(1) mCLK(1) pCLK(1) pcCC(1) pnCC(1) 0*ones(1,2)];

for t = 1:length(mPER)
    tspan = [t t+1]
    Rf = mPERref(t);
    [T,C] = ode45('HU2001_AIFControl_ODE',tspan,C(end,:));
    GeneProteinLevelFull = [GeneProteinLevelFull; C(end,:)];
end

GeneProteinLevelFull = [Cinit; GeneProteinLevelFull(2:end,:)];

%% Plotting

tp = 0:length(mPER)-1;
figure(11)
subplot (4,3,1)
plot(tp,GeneProteinLevelFull(:,1)','LineWidth',2)
hold on
plot(tp,mPERref,'k-','LineWidth',2)
title('mPER')
xlim([0 95])
ylim([0 5])
xticks([0:24:96])
yticks([0:2.5:5])

subplot (4,3,2)
plot(tp,GeneProteinLevelFull(:,2)','LineWidth',2)
title('pPER')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,3)
plot(tp,GeneProteinLevelFull(:,3)','LineWidth',2)
title('mTIM') 
xlim([0 95])
ylim([0 5])
xticks([0:24:96])
yticks([0:2.5:5])

subplot (4,3,4)
plot(tp,GeneProteinLevelFull(:,4)','LineWidth',2)
title('pTIM')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,5)
plot(tp,GeneProteinLevelFull(:,5)','LineWidth',2)
title('pcPT')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])

subplot (4,3,6)
plot(tp,GeneProteinLevelFull(:,6)','LineWidth',2)
title('pnPT')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])

subplot (4,3,7)
plot(tp,GeneProteinLevelFull(:,7)','LineWidth',2)
title('mCLK')
xlim([0 95])
ylim([0 6])
xticks([0:24:96])
yticks([0:3:6])

subplot (4,3,8)
plot(tp,GeneProteinLevelFull(:,8)','LineWidth',2)
title('pCLK')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])


subplot (4,3,9)
plot(tp,GeneProteinLevelFull(:,9)','LineWidth',2)
title('pcCC')
xlim([0 95])
ylim([0 3])
xticks([0:24:96])
yticks([0:1:3])

subplot (4,3,10)
plot(tp,GeneProteinLevelFull(:,10)','LineWidth',2)
title('pnCC')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])

subplot (4,3,11)
plot(tp,GeneProteinLevelFull(:,11)')
title('Z1')
xlim([0 95])
xticks([0:24:96])


subplot (4,3,12)
plot(tp,GeneProteinLevelFull(:,12)')
title('Z2')
xlim([0 95])
xticks([0:24:96])
