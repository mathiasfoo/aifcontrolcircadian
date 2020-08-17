clc
clear all
close all

%% Loading all the experimental WT data under different light condition.


global theta

load drosophila2001LLdata.mat
load drosophila2001referencev1.mat 

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

%% Initial condition

C = [mPER(1) pPER(1) mTIM(1) pTIM(1) pcPT(1) pnPT(1) mCLK(1) pCLK(1) pcCC(1) pnCC(1)];
Cinit = [mPER(1) pPER(1) mTIM(1) pTIM(1) pcPT(1) pnPT(1) mCLK(1) pCLK(1) pcCC(1) pnCC(1)];

for t = 1:length(mPER)
    tspan = [t t+1];
    [T,C] = ode45('HU2001_ESSModel_ODE',tspan,C(end,:));
    GeneProteinLevelFull = [GeneProteinLevelFull; C(end,:)];
end

GeneProteinLevelFull = [Cinit; GeneProteinLevelFull(2:end,:)];

%% Figure Plotting

tp = 0:length(mPER)-1;
figure(1)
subplot (4,3,1)
plot(tp,mPER,tp,GeneProteinLevelFull(:,1)','k--','LineWidth',2)
title('mPER')
xlim([0 95])
ylim([0 5])
xticks([0:24:96])
yticks([0:2.5:5])

subplot (4,3,2)
plot(tp,pPER,tp,GeneProteinLevelFull(:,2)','k--','LineWidth',2)
title('pPER')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,3)
plot(tp,mTIM,tp,GeneProteinLevelFull(:,3)','k--','LineWidth',2)
title('mTIM') 
xlim([0 95])
ylim([0 5])
xticks([0:24:96])
yticks([0:2.5:5])

subplot (4,3,4)
plot(tp,pTIM,tp,GeneProteinLevelFull(:,4)','k--','LineWidth',2)
title('pTIM')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,5)
plot(tp,pcPT,tp,GeneProteinLevelFull(:,5)','k--','LineWidth',2)
title('pcPT')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])

subplot (4,3,6)
plot(tp,pnPT,tp,GeneProteinLevelFull(:,6)','k--','LineWidth',2)
title('pnPT')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])

subplot (4,3,7)
plot(tp,mCLK,tp,GeneProteinLevelFull(:,7)','k--','LineWidth',2)
title('mCLK')
xlim([0 95])
ylim([0 6])
xticks([0:24:96])
yticks([0:3:6])

subplot (4,3,8)
plot(tp,pCLK,tp,GeneProteinLevelFull(:,8)','k--','LineWidth',2)
title('pCLK')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])


subplot (4,3,9)
plot(tp,pcCC,tp,GeneProteinLevelFull(:,9)','k--','LineWidth',2)
title('pcCC')
xlim([0 95])
ylim([0 3])
xticks([0:24:96])
yticks([0:1:3])

subplot (4,3,10)
plot(tp,pnCC,tp,GeneProteinLevelFull(:,10)','k--','LineWidth',2)
title('pnCC')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])
