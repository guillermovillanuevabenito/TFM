format long
close all
clear all
clc

colorA = [0,0,0];
colorB = [1,0.6,0];
colorD = [1,0,1];

%% Load data
load('./data/JYK094_11-RandomReward_03.mat');
addpath('./functions');

lw = 3;
f = 20;

% Range to plot
trange = [-1000 3000]; 
% Range to fit the random lick signal
t0 = -1000; t1 = 3000;
% Licks in trials
licks = data.licks;
% Lick in ITI
l.D = data.licksrandom;
% time 
t = data.lickt;
dt = mean(diff(t));
s.npts = length(t);
rng = t>t0 & t<t1;
% DCN-SNc signal
sig = data.DCN_SNcB;
% DCN-SNc random signal
s.r = data.DCN_SNcBrandom;


% clean up the signals by setting NaNs to 0 and choosing the correct ones
licks = fillmissing(licks,'constant',0); licks = logical(licks);
l.D = fillmissing(l.D,'constant',0); l.D = logical(l.D);
l.DD = l.D(data.logicRandomlick,:);
l.A = licks(data.logicA,:);
l.B = licks(data.logicB,:);
l.D = licks(data.logicD,:);
%
s.r = s.r(data.logicRandomlick,:); s.rm = mean(s.r,1); 
s.A = sig(data.logicA,:); s.Am = mean(s.A,1,'omitnan'); 
s.B = sig(data.logicB,:); s.Bm = mean(s.B,1,'omitnan'); 
s.D = sig(data.logicD,:); s.Dm = mean(s.D,1,'omitnan'); 
s.nr = size(s.r,1);
s.nA = size(s.A,1);
s.nB = size(s.B,1);
s.nD = size(s.D,1);

figure
plot(t,s.Am,'Color',data.colorA); 
hold on
plot(t,s.Bm,'Color',data.colorB); 
plot(t,s.Dm,'Color',data.colorD); 
xlim(trange);  


%{
%% Model Reward
s.Am2fit = s.Am(rng)/(2*max(s.Bm(rng)));
tfit = (0:length(s.Am2fit)-1)*dt;
[s.Akern, s.Agof] = DCN_fit(tfit,s.Am2fit);


t = data.lickt;
figure
%plot(tfit-1000,s.Akern(tfit)*(2*max(s.Am)),'LineWidth',lw,'Color',data.colorA)
hold on
plot(t,s.Am/(2*max(s.Bm(rng))),'Color',data.colorA); 
plot(t,s.Bm/(2*max(s.Bm(rng))),'Color',data.colorB); 
plot(t,s.Dm/(2*max(s.Bm(rng))),'Color',data.colorD); 

xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f)


%%%%% Parameters %%%%%

t_start = 4000;
t_sim = 6000;

t_max = t_start+t_sim;
dt = 0.1;

% Cue
t0_cue = t_start;
t1_cue = t0_cue+250;
J_cue  = 0.08;

% Reward / Higher reward
t0_rew = t_start+150;
t1_rew = t0_rew+250; 
J_rew = 0.06;%0*0.17;
J_hrew = 0.12;%0*0.26;

% Activation function DCN
%m1 = 0.34;
%k1 = 0.15;
%tau1 = 0;

% Activation function recovery variable
%m2 = 0.215;
%k2 = 0.12;
%tau2 = 500;

% Connectivities
R11 = 4;
R12 = -1.8;
R21 = 1;

[DCN_r, d_DCN_r, rec_r, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, J_rew, 0*J_hrew, s.Akern.m1, s.Akern.k1, s.Akern.tau1, s.Akern.m2,s.Akern.k2,s.Akern.tau2, R11, R12, R21);
[DCN_h, d_DCN_h, rec_h, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, 0*J_rew, J_hrew,s.Akern.m1, s.Akern.k1, s.Akern.tau1, s.Akern.m2,s.Akern.k2,s.Akern.tau2, R11, R12, R21);
[DCN_o, d_DCN_o, rec_o, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, -0.06, 0*J_hrew, s.Akern.m1, s.Akern.k1, s.Akern.tau1, s.Akern.m2,s.Akern.k2,s.Akern.tau2, R11, R12, R21);

plot(t-4000,DCN_r,'-','linewidth',2,'Color',data.colorA);
plot(t-4000,DCN_h,'-','linewidth',2,'Color',data.colorB);
plot(t-4000,DCN_o,'-','linewidth',2,'Color',data.colorD);

legend('','','','Reward','Higher-reward','Omission')

%{
plot(t-4000,DCN_r,'-','linewidth',2,'Color',data.colorA);
plot(t-4000,DCN_h,'-','linewidth',2,'Color',data.colorB);
plot(t-4000,DCN_o,'-','linewidth',2,'Color',data.colorD);
legend('','','','Reward','Higher Reward','Omission')
%}


%}
%% Model Higher Reward
s.Bm2fit = s.Bm(rng)/(2*max(s.Bm(rng)));
tfit = (0:length(s.Bm2fit)-1)*dt;
[s.Bkern, s.Bgof] = DCN_fit(tfit,s.Bm2fit);


t = data.lickt;
figure
%plot(tfit-1000,s.Bkern(tfit)*(2*max(s.Bm)),'LineWidth',lw,'Color',data.colorB)
hold on
plot(t,s.Am/(2*max(s.Bm(rng))),'Color',data.colorA); 
plot(t,s.Bm/(2*max(s.Bm(rng))),'Color',data.colorB); 
plot(t,s.Dm/(2*max(s.Bm(rng))),'Color',data.colorD); 

xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f)


%%%%% Parameters %%%%%

t_start = 4000;
t_sim = 6000;

t_max = t_start+t_sim;
dt = 0.1;

% Cue
t0_cue = t_start;
t1_cue = t0_cue+250;
J_cue  = 0.08;

% Reward / Higher reward
t0_rew = t_start+150;
t1_rew = t0_rew+250; 
J_rew = 0.06;%0*0.17;
J_hrew = 0.12;%0*0.26;

% Activation function DCN
%m1 = 0.34;
%k1 = 0.15;
%tau1 = 0;

% Activation function recovery variable
%m2 = 0.215;
%k2 = 0.12;
%tau2 = 500;

% Connectivities
R11 = 4;
R12 = -1.8;
R21 = 1;

[DCN_r, d_DCN_r, rec_r, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, J_rew, 0*J_hrew, s.Bkern.m1, s.Bkern.k1, s.Bkern.tau1, s.Bkern.m2,s.Bkern.k2,s.Bkern.tau2, R11, R12, R21);
[DCN_h, d_DCN_h, rec_h, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, 0*J_rew, J_hrew,s.Bkern.m1, s.Bkern.k1, s.Bkern.tau1, s.Bkern.m2,s.Bkern.k2,s.Bkern.tau2, R11, R12, R21);
[DCN_o, d_DCN_o, rec_o, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, -0.08, 0*J_hrew, s.Bkern.m1, s.Bkern.k1, s.Bkern.tau1, s.Bkern.m2,s.Bkern.k2,s.Bkern.tau2, R11, R12, R21);

plot(t-4000,DCN_r,'-','linewidth',2,'Color',data.colorA);
plot(t-4000,DCN_h,'-','linewidth',2,'Color',data.colorB);
plot(t-4000,DCN_o,'-','linewidth',2,'Color',data.colorD);

legend('','','','Reward','Higher-reward','Omission')




%{
%% Model Omission
s.Dm2fit = s.Dm(rng)/(2*max(s.Dm(rng)));
tfit = (0:length(s.Dm2fit)-1)*dt;
[s.Dkern, s.Dgof] = DCN_fit(tfit,s.Dm2fit);


t = data.lickt;
figure
%plot(tfit-1000,s.Dkern(tfit)*(2*max(s.Dm)),'LineWidth',lw,'Color',data.colorD)
hold on
plot(t,s.Am,'Color',data.colorA); 
plot(t,s.Bm,'Color',data.colorB); 
plot(t,s.Dm,'Color',data.colorD); 

xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f)


%%%%% Parameters %%%%%

t_start = 4000;
t_sim = 6000;

t_max = t_start+t_sim;
dt = 0.1;

% Cue
t0_cue = t_start;
t1_cue = t0_cue+250;
J_cue  = 0.2;

% Reward / Higher reward
t0_rew = t_start+150;
t1_rew = t0_rew+250; 
J_rew = 0.18;
J_hrew = 0.2;

% Activation function DCN
%m1 = 0.34;
%k1 = 0.15;
%tau1 = 0;

% Activation function recovery variable
%m2 = 0.215;
%k2 = 0.12;
%tau2 = 500;

% Connectivities
R11 = 4;
R12 = -1.8;
R21 = 1;

[DCN_r, d_DCN_r, rec_r, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, J_rew, 0*J_hrew, s.Dkern.m1, s.Dkern.k1, s.Dkern.tau1, s.Dkern.m2,s.Dkern.k2,s.Dkern.tau2, R11, R12, R21);
[DCN_h, d_DCN_h, rec_h, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, 0*J_rew, J_hrew,s.Dkern.m1, s.Dkern.k1, s.Dkern.tau1, s.Dkern.m2,s.Dkern.k2,s.Dkern.tau2, R11, R12, R21);
[DCN_o, d_DCN_o, rec_o, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, 0*J_rew, 0*J_hrew, s.Dkern.m1, s.Dkern.k1, s.Dkern.tau1, s.Dkern.m2,s.Dkern.k2,s.Dkern.tau2, R11, R12, R21);

plot(t-4000,DCN_r*(2*max(s.Dm)),'-','linewidth',2,'Color',data.colorA);
plot(t-4000,DCN_h*(2*max(s.Dm)),'-','linewidth',2,'Color',data.colorB);
plot(t-4000,DCN_o*(2*max(s.Dm)),'-','linewidth',2,'Color',data.colorD);
%}