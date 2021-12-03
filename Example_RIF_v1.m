%  Test Example on Bat sound by Patrick Flandrin
% 
%  Ref: G. Barbarino and A. Cicone. 
%  "Stabilization and Variations to the Adaptive Local Iterative Filtering Algorithm: the Fast Resampled Iterative Filtering Method". 
%  Submitted 2021
%  arXiv: http://arxiv.org/abs/2111.02764

set(0,'defaultTextInterpreter','latex');
load('batrich.mat')
figure
fs=60;
plot((1:length(x))/fs,x)

axis tight
set(gcf,'units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',36,'DefaultLineLineWidth',2)
set(gca,'FontSize',36)

%% FIF decomposition
opts=Settings_FIF_v3('alpha',10,'NIMFs',19);
[IMFs,stats] = FIF_v2_12(x,opts);

fig=plot_imf_v10(IMFs,(1:length(x))/fs,20)
set(fig,'units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',36,'DefaultLineLineWidth',2)


%% IMFogram time-frequency representation

IMFogram_v1(IMFs,fs,90,[],100,15);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

%% Inst Freqs - 1st curve

t11 = 2.9*60; 
y11 = 27;
t12 = 6.1*60;
y12 = 14;

t = 1:400; % fs=60 Hz
if1 = (y12-y11)/(t12-t11)*(t-t11)+y11;

%% Inst Freqs - 2nd curve 
% Using cftool for curve fitting. We choose exponential
t2 = [0.6 1.9 2.9 3.8 4.8 5.8]*60; 
y2 = [28.61 22.6 18.7 15.7 12.8 10.3];
% General model Exp1:
%     f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
       a =       32.32;%  (31.57, 33.07)
       b =   -0.003203;%  (-0.003345, -0.00306)
% Goodness of fit:
%   SSE: 0.2099
%   R-square: 0.9991
%   Adjusted R-square: 0.9988
%   RMSE: 0.2291
if2 = a*exp(b*t);
%% Inst Freqs - 3rd curve 
% Using cftool for curve fitting. We choose exponential
t3 = [0.06 0.38 0.833 1.25  1.733 2.27 2.583 3.05 3.533 4.267 4.917 5.567 6.4]*60; 
y3 = [21.7 15.8 13.31 11.94 10.76 9.81 9.309 8.776 7.858 6.822 6.142 5.49 5];
%General model Exp2:
%     f(x) = a*exp(b*x) + c*exp(d*x)
%Coefficients (with 95% confidence bounds):
       a =       8.909;%  (8.346, 9.473)
       b =    -0.06368;%  (-0.07411, -0.05325)
       c =       14.76;%  (14.37, 15.15)
       d =   -0.002941;%  (-0.003083, -0.002798)

% Goodness of fit:
%   SSE: 0.1738
%   R-square: 0.9993
%   Adjusted R-square: 0.9991
%   RMSE: 0.139

if3 = a*exp(b*t) + c*exp(d*t);
%% Inst freq curves

figure('units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',36,'DefaultLineLineWidth',2)
t=(1:400)/60;
plot(t,if1)
hold on
plot(t,if2)
plot(t,if3)
axis([1/60 386/60 0 30])

%%
dt=1/60;
L1 = 1./if1/dt;
L2 = 1./if2/dt;
L3 = 1./if3/dt;
% figure
% plot(L1)
% hold on
% plot(L2)
% plot(L3)

%% FRIF method after reallignement of the peaks in Fourier domain

opts = Settings_FRIF_v1('MaxInner',1000,'Xi',1.5,'delta',10^-6,'alpha',30,'plots',0,'UpSampling',128);
[IMFs2,stats] = FRIF_v1_3(x,[L1],opts,2);

opts = Settings_FRIF_v1('MaxInner',1000,'Xi',2.1,'delta',10^-6,'alpha',30,'plots',0,'UpSampling',128);
[IMFs3,stats] = FRIF_v1_3(IMFs2(end,:),[L2],opts,2);

opts = Settings_FRIF_v1('MaxInner',1000,'Xi',2.25,'delta',10^-6,'alpha',30,'plots',0,'UpSampling',128);
[IMFs4,stats] = FRIF_v1_3(IMFs3(end,:),[L3],opts,4);

IMFs_new=[IMFs2(1,:);IMFs3(1,:);IMFs4];
plot_imf_v10(IMFs_new,(1:length(x))/fs,4)
set(gcf,'units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',36,'DefaultLineLineWidth',2)

%% TFR

[X,Y,A,IF,fig]=IMFogram_v1(IMFs_new(1,:),fs,90,[],100,15);

[X,Y,A,IF,fig]=IMFogram_v1(IMFs_new(2,:),fs,90,[],100,15);

[X,Y,A,IF,fig]=IMFogram_v1(IMFs_new(3,:),fs,90,[],100,15);


