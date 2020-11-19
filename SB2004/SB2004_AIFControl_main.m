clc
clear all
close all

%% Loading all the experimental WT data under different light condition.

global theta Rf eta theta1 theta2 gammaC

eta = 15.2; 
theta1 = 3;
theta2 = 1;
gammaC = 1.5;

load simplemammalianLLdata.mat
load mammalianreference.mat


%% Model Parameters
prmPER2CRY = [0.14959,-5.5475,0.045101,0.15002];
prpcPER2CRY = [0.32152,1.7465,0.27062];
prpnPER2CRY = [0.17611,1.3147,0.10316];
prmBMAL1 = [0.80972,2.057,1.5816];
prpcBMAL1 = [0.20242,0.95272,0.34725];
prpnBMAL1 = [0.38498,1.0573,0.22119];
prpaBMAL1 = [0.11685,1.0321,0.12062];

nonlineartheta = [prmPER2CRY prpcPER2CRY prpnPER2CRY prmBMAL1 prpcBMAL1 prpnBMAL1 prpaBMAL1];

theta = nonlineartheta;
GeneProteinLevelFull = [];

%% Initial condition

C = [mPER2CRY(1) pcPER2CRY(1) pnPER2CRY(1) mBMAL1(1) pcBMAL1(1) pnBMAL1(1) paBMAL1(1) 0*ones(1,2)];
Cinit = [mPER2CRY(1) pcPER2CRY(1) pnPER2CRY(1) mBMAL1(1) pcBMAL1(1) pnBMAL1(1) paBMAL1(1) 0*ones(1,2)];

for t = 1:length(mPER2CRY)
    tspan = [t t+1]
    Rf = mBMAL1ref(t);
    [T,C] = ode45('SB2004_AIFControl_ODE',tspan,C(end,:));
    GeneProteinLevelFull = [GeneProteinLevelFull; C(end,:)];
end

GeneProteinLevelFull = [Cinit; GeneProteinLevelFull(2:end,:)];

%% Figure Plotting

tp = 0:length(mPER2CRY)-1;
figure(1)
subplot (4,3,1)
plot(tp,GeneProteinLevelFull(:,1)','LineWidth',2)
title('mPER2CRY')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,2)
plot(tp,GeneProteinLevelFull(:,2)','LineWidth',2)
title('pcPER2CRY')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,3)
plot(tp,GeneProteinLevelFull(:,3)','LineWidth',2)
title('pnPER2CRY') 
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,4)
plot(tp,GeneProteinLevelFull(:,4)','LineWidth',2)
hold on
plot(tp,mBMAL1ref,'k-','LineWidth',2)
title('mBMAL1')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,5)
plot(tp,GeneProteinLevelFull(:,5)','LineWidth',2)
title('pcBMAL1')
xlim([0 95])
ylim([0 1])
xticks([0:24:96])
yticks([0:0.5:1])

subplot (4,3,6)
plot(tp,GeneProteinLevelFull(:,6)','LineWidth',2)
title('pnBMAL1')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,7)
plot(tp,GeneProteinLevelFull(:,7)','LineWidth',2)
title('paBMAL1')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,8)
plot(tp,GeneProteinLevelFull(:,8)','LineWidth',2)
title('Z1')
xlim([0 95])
xticks([0:24:96])

subplot (4,3,9)
plot(tp,GeneProteinLevelFull(:,9)','LineWidth',2)
title('Z2')
xlim([0 95])
xticks([0:24:96])
