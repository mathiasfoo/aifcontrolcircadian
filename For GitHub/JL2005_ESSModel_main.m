clc
clear all
close all

%% Loading training data from JL2005 using Michaelis Menten model


global EL theta

load LockeJTB_training.mat

r1 = 1; r2 = 96;
mLHY = mLHY12(r1:r2);
pcLHY = pcLHY12(r1:r2);
pnLHY = pnLHY12(r1:r2);
pLHY = pcLHY + pnLHY;
mTOC1 = mTOC112(r1:r2);
pcTOC1 = pcTOC112(r1:r2);
pnTOC1 = pnTOC112(r1:r2);
pP = pP12(r1:r2);
TL = TL12(r1:r2);

% Model Parameters
prmLHY = [4.7463,1.4361,1.3276,1.2875];
prpcLHY = [0.4074,2.0764,1.7443];
prpnLHY = [0.9147,1.7344,1.4088];
prmTOC1 = [0.0484,-1.4874,1.5163];
prpcTOC1 = [2.3917,1.9354,0.2075];
prpnTOC1 = [0.0378,1.1309,0.0268];
prpP = [0.4080,0.4886,1.5690];

nonlineartheta = [prmLHY prpcLHY prpnLHY prmTOC1 prpcTOC1 prpnTOC1 prpP];

theta = nonlineartheta;

GeneProteinLevelFull = [];

% return

%% Initial condition

C = [mLHY(1) pcLHY(1) pnLHY(1) mTOC1(1) pcTOC1(1) pnTOC1(1) pP(1)];
Cinit = [mLHY(1) pcLHY(1) pnLHY(1) mTOC1(1) pcTOC1(1) pnTOC1(1) pP(1)];

for t = 1:length(TL)
    tspan = [t t+1];
    EL = TL(t);    
    [T,C] = ode45('JL2005_ESSModel_ODE',tspan,C(end,:));
    GeneProteinLevelFull = [GeneProteinLevelFull; C(end,:)];
end

GeneProteinLevelFull = [Cinit; GeneProteinLevelFull(1:end-1,:)];

%% Plotting

tp = 0:length(TL)-1;

figure(1)
subplot(3,3,1)
plot(tp,mLHY,tp,GeneProteinLevelFull(:,1),'k--','LineWidth',2)
title('mLHY')
xlim([0 95])
ylim([0 3])
xticks([0:24:96])
yticks([0:1:3])

subplot(3,3,2)
plot(tp,pcLHY,tp,GeneProteinLevelFull(:,2)','k--','LineWidth',2)
title('pcLHY')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot(3,3,3)
plot(tp,pnLHY,tp,GeneProteinLevelFull(:,3)','k--','LineWidth',2)
title('pnLHY') 
xlim([0 95])
ylim([0 1])
xticks([0:24:96])
yticks([0:0.5:1])

subplot(3,3,4)
plot(tp,mTOC1,tp,GeneProteinLevelFull(:,4)','k--','LineWidth',2)
title('mTOC1')
xlim([0 95])
ylim([0 0.5])
xticks([0:24:96])
yticks([0:0.25:0.5])

subplot(3,3,5)
plot(tp,pcTOC1,tp,GeneProteinLevelFull(:,5)','k--','LineWidth',2)
title('pcTOC1')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot(3,3,6)
plot(tp,pnTOC1,tp,GeneProteinLevelFull(:,6)','k--','LineWidth',2)
title('pnTOC1')
xlim([0 95])
ylim([1 1.5])
xticks([0:24:96])
yticks([1:0.25:1.5])

subplot(3,3,7)
plot(tp,pP,tp,GeneProteinLevelFull(:,7)','k--','LineWidth',2)
title('pP')
xlim([0 95])
ylim([0 1])
xticks([0:24:96])
yticks([0:0.5:1])
