clc
clear all
close all

% De Caluwe et al (2016) FPS

global theta EL



load DeCaluwe2016FPS12L12DLLdata.mat
load JD2016referencesim.mat

r1 = 1; r2 = 96;
mLHY = mLHY12(r1:r2)';
pLHY = pLHY12(r1:r2)';
mP97 = mP9712(r1:r2)';
pP97 = pP9712(r1:r2)';
mP51 = mP5112(r1:r2)';
pP51 = pP5112(r1:r2)';
mEL = mEL12(r1:r2)';
pEL = pEL12(r1:r2)';
pP = pP12(r1:r2)';
mPIF = mPIF12(r1:r2)';
pPIF = pPIF12(r1:r2)';
pHYP = phHYP12(r1:r2)';
TL = TL12(r1:r2);

%% Model Parameters
prmLHY = [0.72294,-1.2502,-0.042071,0.75237,-0.97172];
prpLHY = [0.80635,0.74888,0.74688];
prmP97 = [0.026394,0.12407,-2.3687,-0.10561,0.35392,1.5597];
prpP97 = [1.0733,0.29692,-0.050611];
prmP51 = [0.38365,-0.71673,0.17357,1.013];
prpP51 = [0.55202,0.41674,-0.27927];
prmEL = [0.010977,-1.8819,-1.3932,-0.017199,0.48139];
prpEL = [1.0696,1.3282,0.91055];
prpP = [0.34819,0.36754,-0.64076];
prmPIF = [0.036342,-0.28081,0.19706];
prpPIF = [0.30969,0.18321,-2.0456];
prpHYP = [0.17399,1.1843];

linearthetanl = [prmLHY prpLHY prmP97 prpP97 prmP51 prpP51 prmEL prpEL prpP prmPIF prpPIF prpHYP];

theta = linearthetanl;

C = [mLHY(1) pLHY(1) mP97(1) pP97(1) mP51(1) pP51(1) mEL(1) pEL(1) pP(1) mPIF(1) pPIF(1) pHYP(1)];
Cinit = [mLHY(1) pLHY(1) mP97(1) pP97(1) mP51(1) pP51(1) mEL(1) pEL(1) pP(1) mPIF(1) pPIF(1) pHYP(1)];

ProteinLevel = [];

for t = 1:length(TL)
    tspan = [t t+1];
    EL = TL(t);    
    [T,C] = ode45('JD2016_ESSModel_ODE',tspan,C(end,:));
    ProteinLevel = [ProteinLevel; C(end,:)];
end
ProteinLevel = [Cinit; ProteinLevel(2:end,:)];

%% Figure plotting
tp = 0:length(TL)-1;

figure(1)
subplot(4,4,1)
plot(tp,mLHY,tp,ProteinLevel(:,1),'k--','LineWidth',2)
title('mLHY')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot(4,4,2)
plot(tp,pLHY,tp,ProteinLevel(:,2),'k--','LineWidth',2)
title('pcLHY')
xlim([0 95])
ylim([0 3])
xticks([0:24:96])
yticks([0:1:3])

subplot(4,4,3)
plot(tp,mP97,tp,ProteinLevel(:,3),'k--','LineWidth',2)
title('mP97')
xlim([0 95])
ylim([0 4])
xticks([0:24:96])
yticks([0:2:4])

subplot(4,4,4)
plot(tp,pP97,tp,ProteinLevel(:,4),'k--','LineWidth',2)
title('pP97')
xlim([0 95])
ylim([0 10])
xticks([0:24:96])
yticks([0:5:10])

subplot(4,4,5)
plot(tp,mP51,tp,ProteinLevel(:,5),'k--','LineWidth',2)
title('mP51')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot(4,4,6)
plot(tp,pP51,tp,ProteinLevel(:,6),'k--','LineWidth',2)
title('pP51')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot(4,4,7)
plot(tp,mEL,tp,ProteinLevel(:,7),'k--','LineWidth',2)
title('mEL')
xlim([0 95])
ylim([0 2])
xticks([0:24:96])
yticks([0:1:2])

subplot(4,4,8)
plot(tp,pEL,tp,ProteinLevel(:,8),'k--','LineWidth',2)
title('pEL')
xlim([0 95])
ylim([0 3])
xticks([0:24:96])
yticks([0:1:3])
 
subplot(4,4,9)
plot(tp,pP,tp,ProteinLevel(:,9),'k--','LineWidth',2)
title('P')
xlim([0 95])
ylim([0 1])
xticks([0:24:96])
yticks([0:0.5:1])

subplot(4,4,10)
plot(tp,mPIF,tp,ProteinLevel(:,10),'k--','LineWidth',2)
title('mPIF')
xlim([0 95])
ylim([0 1])
xticks([0:24:96])
yticks([0:0.5:1])

subplot(4,4,11)
plot(tp,pPIF,tp,ProteinLevel(:,11),'k--','LineWidth',2)
title('pPIF')
xlim([0 95])
ylim([0 1])
xticks([0:24:96])
yticks([0:0.5:1])

subplot(4,4,12)
plot(tp,ProteinLevel(:,12),'LineWidth',2)
hold on
plot(tp,pHYP,'k-','LineWidth',2)
hold on
title('pHYP')
xlim([0 95])
ylim([6 9])
xticks([0:24:96])
yticks([6:1:9])
