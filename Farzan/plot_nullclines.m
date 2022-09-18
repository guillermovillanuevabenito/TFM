format long
close all
clear all
clc

lw = 3;
f = 20;
load('./data/JYK094_11-RandomReward_03.mat');

%% Load DCN model

%%%%% Parameters %%%%%

% Range to plot
trange = [-1000 3000]; 

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

[DCN_r, d_DCN_r, rec_r, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, J_rew, 0*J_hrew, m1, k1, tau1, m2,k2,tau2, R11, R12, R21);
[DCN_h, d_DCN_h, rec_h, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, 0*J_rew, J_hrew,m1, k1, tau1, m2,k2,tau2, R11, R12, R21);
[DCN_o, d_DCN_o, rec_o, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, -J_rew, 0*J_hrew, m1, k1, tau1, m2,k2,tau2,R11, R12, R21);

plot(t-4000,DCN_r,'-','linewidth',lw,'Color',data.colorA);
plot(t-4000,DCN_h,'-','linewidth',lw,'Color',data.colorB);
plot(t-4000,DCN_o,'-','linewidth',lw,'Color',data.colorD);
xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f)
legend('Reward','Higher-reward','Omission')

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


plot(DCN_r(t>=4000),rec_r(t>=4000),'-','linewidth',2,'Color',data.colorA);
plot(DCN_o(t>=4000),rec_o(t>=4000),'-','linewidth',2,'Color',data.colorD);
plot(DCN_h(t>=4000),rec_h(t>=4000),'-','linewidth',2,'Color',data.colorB);

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
