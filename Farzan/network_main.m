format long
close all
clear all
clc
%%

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


t = data.lickt;
figure
hold on

plot(t,s.Am/(2*max(s.Bm(rng))),'Color',data.colorA); 
plot(t,s.Bm/(2*max(s.Bm(rng))),'Color',data.colorB); 
plot(t,s.Dm/(2*max(s.Bm(rng))),'Color',data.colorD); 


xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f)

%% Load DCN model

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
J_rew = 0.08; %0*0.17;
J_hrew = 0.15; %0*0.26;

% Activation function DCN
m1 = 0.484;
k1 = 0.243;
tau1 = 433;

% Activation function recovery variable
m2 = 0.293;
k2 = 0.1476;
tau2 = 108;

% Connectivities
R11 = 4;
R12 = -1.8;
R21 = 1;

figure
hold on
plot(t,s.Am/(2*max(s.Bm(rng))),'Color',data.colorA); 
plot(t,s.Bm/(2*max(s.Bm(rng))),'Color',data.colorB); 
plot(t,s.Dm/(2*max(s.Bm(rng))),'Color',data.colorD); 


[DCN_r, d_DCN_r, rec_r, tt] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, J_rew, 0*J_hrew, m1, k1, tau1, m2,k2,tau2, R11, R12, R21);
[DCN_h, d_DCN_h, rec_h, tt] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, 0*J_rew, J_hrew,m1, k1, tau1, m2,k2,tau2, R11, R12, R21);
[DCN_o, d_DCN_o, rec_o, tt] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, -J_rew, 0*J_hrew, m1, k1, tau1, m2,k2,tau2,R11, R12, R21);

plot(tt-4000,DCN_r,'-','linewidth',lw,'Color',data.colorA);
plot(tt-4000,DCN_h,'-','linewidth',lw,'Color',data.colorB);
plot(tt-4000,DCN_o,'-','linewidth',lw,'Color',data.colorD);
xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('','','','Reward','Higher-reward','Omission')


figure
hold on
plot(tt-4000,rec_r,'-','linewidth',lw,'Color',data.colorA);
plot(tt-4000,rec_h,'-','linewidth',lw,'Color',data.colorB);
plot(tt-4000,rec_o,'-','linewidth',lw,'Color',data.colorD);
xlim(trange);  
ylabel('a');
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Reward','Higher-reward','Omission')

%% Plot sensorimotor inputs

pulse = @(t,t0,t1) heaviside(t-t0).*heaviside(t1-t);
figure

subplot(3,1,1);
hold on
plot(tt-4000,J_cue*pulse(tt,t0_cue,t1_cue)-J_rew*pulse(tt,t0_rew,t1_rew),'-','linewidth',lw,'Color',data.colorD);

xlim(trange);  
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Omission')

subplot(3,1,2);
hold on
plot(tt-4000,J_cue*pulse(tt,t0_cue,t1_cue) + J_rew*pulse(tt,t0_rew,t1_rew),'-','linewidth',lw,'Color',data.colorA);

xlim(trange);  
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Reward')

subplot(3,1,3);
hold on
plot(tt-4000,J_cue*pulse(tt,t0_cue,t1_cue)+J_hrew*pulse(tt,t0_rew,t1_rew),'-','linewidth',lw,'Color',data.colorB);

xlim(trange);  
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Higher-reward')


figure
subplot(2,1,1);
hold on
plot(tt-4000,J_cue*pulse(tt,t0_cue,t1_cue),'-','linewidth',lw,'Color','blue');
xlim(trange);  
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Cue')

subplot(2,1,2);
hold on
plot(tt-4000,J_hrew*pulse(tt,t0_rew,t1_rew),'-','linewidth',lw,'Color',data.colorB);
plot(tt-4000,-J_rew*pulse(tt,t0_rew,t1_rew),'-','linewidth',lw,'Color',data.colorD);
plot(tt-4000,J_rew*pulse(tt,t0_rew,t1_rew),'-','linewidth',lw,'Color',data.colorA);

xlim(trange);  
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Higher-reward','Omission','Reward')


%% Plot nullclines
figure
subplot(1,2,1);
I = 0.0;

x_sup = (2*k1+m1-I-R11)/R12;
x_inf = (m1-I-2*k1)/R12;


xx2 = -1:0.0001:x_sup;
xx = x_inf:0.0001:x_sup;
yy = w_nullcline(xx2,m1,k1,m2,k2,R11,R12,R21);
plot(xx2,yy,'linewidth',1,'Color','r')
hold on
yy = z_nullcline(xx,m1,k1,m2,k2,R11,R12,R21,I);
plot(yy,xx,'linewidth',1,'Color','b')
plot([1 1],[x_sup -1],'linewidth',1,'Color','b')
plot([0 0],[x_inf 2],'linewidth',1,'Color','b')


plot(DCN_r(t>=4000),rec_r(t>=4000),'-','linewidth',2,'Color',colorA);
plot(DCN_o(t>=4000),rec_o(t>=4000),'-','linewidth',2,'Color',colorD);
plot(DCN_h(t>=4000),rec_h(t>=4000),'-','linewidth',2,'Color',colorB);

plot(1,1,'.k','MarkerSize',20,'linewidth',2)
plot(0.593,1,'ok','MarkerSize',6,'linewidth',2)
plot(0,0.004,'.k','MarkerSize',20,'linewidth',2)
legend('','','','','Reward trajectory','Omission trajectory','Higher-reward trajectory','','','','Location','northoutside')

%legend('a-nullcline','r_{1e}-nullcline','','','Location','northoutside')
%title('Phase plane')
set(gca,'fontsize',24);
xlabel('r_{1e}');
ylabel('a');
xlim([0-0.1 1+0.1])
ylim([0-0.1 1+0.1])


%{
a_t = -0.5:0.1:1.5;
z_t = -0.5:0.1:1.5;

for ii=1:length(a_t)
    for jj=1:length(z_t)
        [DCN_r, d_DCN_r, rec_r, t] = rate_model_2d_test(t_max/3, dt, t0_cue, t1_cue, 0*J_cue,t0_rew, t1_rew, 0*J_rew, 0*J_hrew, m1, k1, tau1, m2,k2,tau2, R11, R12, R21,z_t(jj),a_t(ii));
        plot(DCN_r,rec_r,'-','linewidth',1)
    end
end

%}



subplot(1,2,2);
I = 0.15;

x_sup = (2*k1+m1-I-R11)/R12;
x_inf = (m1-I-2*k1)/R12;


xx2 = -1:0.0001:x_sup;
xx = x_inf:0.0001:x_sup;
yy = w_nullcline(xx2,m1,k1,m2,k2,R11,R12,R21);
plot(xx2,yy,'linewidth',1,'Color','r')
hold on
yy = z_nullcline(xx,m1,k1,m2,k2,R11,R12,R21,I);
plot(yy,xx,'linewidth',1,'Color','b')
plot([1 1],[x_sup -1],'linewidth',1,'Color','b')
plot([0 0],[x_inf 2],'linewidth',1,'Color','b')
plot(1,1,'.k','MarkerSize',20,'linewidth',2)

legend('a-nullcline','r_{1e}-nullcline','','','','Location','northoutside')
%title('Phase plane')
set(gca,'fontsize',24);
xlabel('r_{1e}');
ylabel('a');
xlim([0-0.1 1+0.1])
ylim([0-0.1 1+0.1])

%% DCN-SNc Network

colorA = [0,0,0];
colorB = [1,0.6,0];
colorD = [1,0,1];

% Areas:
%       1: DCN
%       2: SNc
%       3: SNr
%
% Neurons (populations, nodes)
%
% DCN: x1e
% SNc: x2e, x2i
% SNr: x3i
% Motor: x4e

% [t] = ms

% Connectivity convention: Jxy: J:y->x

% Parameters

x2ethr = -0.2;
x2eslp = 1.4;

x2ithr = -0.01;
x2islp = 1.45;

x3ithr = -0.06;
x3islp = 0.3;

tau2e = 410; %250
tau2i = 220; % 100
tau3i = 120; %170


J2i1e = 0;

J2e1e = 0.89;    % 0.9
J3i1e = 0.21;  % 0.5

J2i2e = 0.1*2.7;  %1
J2e2i = 0.2*4.3;  %1
J2e3i = 0.2*1.001;  % 0.63

J1e4e = 0;
J4e1e = 0;


t_start = 4000;
t_sim = 6000;

t_max = t_start+t_sim;

Tmax = t_max;
dt = 0.1;
t = 0:dt:Tmax;

% Inputs

%2e

ton2ecue = t_start;
toff2ecue = ton2ecue + 250;
A2ecue = 0.04;

load('./data/JYK0942_11-RandomReward_03.mat');
% Range to fit the random lick signal
t0 = -1000; t1 = 3000;
% Licks in trials
licks = data.licks;

% time 
tt = data.lickt;
dtt = mean(diff(tt));
s.npts = length(tt);
rng = tt>t0 & tt<t1;
% DCN-SNc signal
sig = data.SNcA;

% clean up the signals by setting NaNs to 0 and choosing the correct ones
licks = fillmissing(licks,'constant',0); licks = logical(licks);
l.A = licks(data.logicA,:);
l.B = licks(data.logicB,:);
l.D = licks(data.logicD,:);
%
s.A = sig(data.logicA,:); s.Am = mean(s.A,1,'omitnan'); 
s.B = sig(data.logicB,:); s.Bm = mean(s.B,1,'omitnan'); 
s.D = sig(data.logicD,:); s.Dm = mean(s.D,1,'omitnan'); 
s.nA = size(s.A,1);
s.nB = size(s.B,1);
s.nD = size(s.D,1);


s.Am2fit = s.Am(rng)/(2*max(s.Bm(rng)));
tfit = (0:length(s.Am2fit)-1)*dtt;


[SNC.DA_o,x2i_o, x3i_o] = rate_model_network(s.Am2fit,DCN_o,Tmax,dt,x2ethr,x2eslp,x2ithr,x2islp,x3ithr,x3islp,tau2e,tau2i,tau3i,J2i1e,J2e1e,J2i2e,J2e2i,J2e3i,J3i1e,ton2ecue,toff2ecue,A2ecue);
[SNC.DA_r,x2i_r, x3i_r] = rate_model_network(s.Am2fit,DCN_r,Tmax,dt,x2ethr,x2eslp,x2ithr,x2islp,x3ithr,x3islp,tau2e,tau2i,tau3i,J2i1e,J2e1e,J2i2e,J2e2i,J2e3i,J3i1e,ton2ecue,toff2ecue,A2ecue);
[SNC.DA_h,x2i_h, x3i_h] = rate_model_network(s.Am2fit,DCN_h,Tmax,dt,x2ethr,x2eslp,x2ithr,x2islp,x3ithr,x3islp,tau2e,tau2i,tau3i,J2i1e,J2e1e,J2i2e,J2e2i,J2e3i,J3i1e,ton2ecue,toff2ecue,A2ecue);


figure
hold on
plot(tfit-1000,SNC.DA_r-SNC.DA_r(17),'LineWidth',lw,'Color',colorA)
plot(tfit-1000,SNC.DA_h-SNC.DA_r(17),'LineWidth',lw,'Color',colorB)
plot(tfit-1000,SNC.DA_o-SNC.DA_r(17),'LineWidth',lw,'Color',colorD)
kk= 70;

plot(tt,s.Am/(kk*max(s.Bm)),'Color',data.colorA); 
hold on;
plot(tt,s.Bm/(kk*max(s.Bm)),'Color',data.colorB); 
plot(tt,s.Dm/(kk*max(s.Bm)),'Color',data.colorD); 
legend('Reward','Higher-reward','Omission','','','')


set(gca,'fontsize',24);
%axis([0 Tmax 0 0.5]);
xlim([-1000 3000])
xlabel('t[ms]');
ylabel('SNc.DA');




figure
plot(tfit-1000,x2i_r,'LineWidth',lw,'Color',colorA)
hold on
plot(tfit-1000,x2i_h,'LineWidth',lw,'Color',colorB)
plot(tfit-1000,x2i_o,'LineWidth',lw,'Color',colorD)
legend('Reward','Higher-reward','Omission')

set(gca,'fontsize',24);

xlim([-1000 3000])
xlabel('t[ms]');
ylabel('SNc.GABA');

figure
plot(tfit-1000,x3i_r,'LineWidth',lw,'Color',colorA)
hold on
plot(tfit-1000,x3i_h,'LineWidth',lw,'Color',colorB)
plot(tfit-1000,x3i_o,'LineWidth',lw,'Color',colorD)
legend('Reward','Higher-reward','Omission')
xlim([-1000 3000])
set(gca,'fontsize',24);

xlabel('t[ms]');
ylabel('SNr.GABA');