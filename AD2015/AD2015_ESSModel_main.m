clc
clear all
close all

%% Loading all the experimental WT data under different light condition.

global theta

load neurosporaLLdata.mat
load neurosporareferenceSS.mat

%% Model Parameters
prmFRQ = [0.80569,3.9937,0.27491];
prpcFRQ = [1.5144,0.25105];
prpnFRQ = [0.050751,0.28667,48.8696,1.0255];
prmWC1 = [0.43796,-0.12732,0.11805];
prpcWC1 = [0.066412,1.0938,0.20418,64.3947];
prpnWC1 = [42.0093,0.014821,51.0546,0.90579];
prpFRQWC1 = [52.0211,6.3082];
prmCSP1 = [0.063,0.92623,-0.17928,1.0502];
prpCSP1 = [52.412,1.4815];

nonlineartheta = [prmFRQ prpcFRQ prpnFRQ prmWC1 prpcWC1 prpnWC1 prpFRQWC1 prmCSP1 prpCSP1];

theta = nonlineartheta;
GeneProteinLevelFull = [];

%% Initial condition

C = [mFRQ(1) pcFRQ(1) pnFRQ(1) mWC1(1) pcWC1(1) pnWC1(1) pFRQWC1(1) mCSP1(1) pCSP1(1)];
Cinit = [mFRQ(1) pcFRQ(1) pnFRQ(1) mWC1(1) pcWC1(1) pnWC1(1) pFRQWC1(1) mCSP1(1) pCSP1(1)];

for t = 1:length(mFRQ)
    tspan = [t t+1]
    [T,C] = ode45('AD2015_ESSModel_ODE',tspan,C(end,:));
    GeneProteinLevelFull = [GeneProteinLevelFull; C(end,:)];
end

GeneProteinLevelFull = [Cinit; GeneProteinLevelFull(2:end,:)];

%% Figure Plotting

tp = 0:length(mFRQ)-1;

figure(11)
subplot (4,3,1)
plot(tp,mFRQ,tp,GeneProteinLevelFull(:,1)','k--','LineWidth',2)
hold on
title('mFRQ')
xlim([0 95])
ylim([0 3])
xticks([0:24:96])
yticks([0:1:3])

subplot (4,3,2)
plot(tp,pcFRQ,tp,GeneProteinLevelFull(:,2)','k--','LineWidth',2)
hold on
title('pcFRQ')
xlim([0 95])
ylim([0 12])
xticks([0:24:96])
yticks([0:6:12])

subplot (4,3,3)
plot(tp,pnFRQ,tp,GeneProteinLevelFull(:,3)','k--','LineWidth',2)
hold on
title('pnFRQ') 
xlim([0 95])
ylim([0 0.5])
xticks([0:24:96])
yticks([0:0.25:0.5])

subplot (4,3,4)
plot(tp,mWC1,tp,GeneProteinLevelFull(:,4)','k--','LineWidth',2)
hold on
title('mWC1')
xlim([0 95])
ylim([3 5])
xticks([0:24:96])
yticks([3:1:5])

subplot (4,3,5)
plot(tp,pcWC1,tp,GeneProteinLevelFull(:,5)','k--','LineWidth',2)
hold on
title('pcWC1')
xlim([0 95])
ylim([4e-3 8e-3])
xticks([0:24:96])
yticks([4e-3:2e-3:8e-3])

subplot (4,3,6)
plot(tp,pnWC1,tp,GeneProteinLevelFull(:,6)','k--','LineWidth',2)
hold on
title('pnWC1')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot (4,3,7)
plot(tp,pFRQWC1,tp,GeneProteinLevelFull(:,7)','k--','LineWidth',2)
hold on
title('pFRQWC1')
xlim([0 95])
ylim([0 0.12])
xticks([0:24:96])
yticks([0:0.06:0.12])

subplot (4,3,8)
plot(tp,mCSP1,tp,GeneProteinLevelFull(:,8)','k--','LineWidth',2)
hold on
title('mCSP1')
xlim([0 95])
ylim([0 0.06])
xticks([0:24:96])
yticks([0:0.03:0.06])

subplot (4,3,9)
plot(tp,pCSP1,tp,GeneProteinLevelFull(:,9)','k--','LineWidth',2)
hold on
title('pCSP1')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])
