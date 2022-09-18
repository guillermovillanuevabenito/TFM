format long
close all
clearvars

%% Load data
load('./data/JYK094_11-RandomReward_03.mat');
addpath('./functions');

lw = 3;
f = 20;

% Range to plot
trange = [-1000 3000]; 
% Range to fit the random lick signal
t0 = -250; t1 = 500;
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

% Average plot
figure
plot(t,s.Am,'LineWidth',lw,'Color',data.colorA); 
hold on;
plot(t,s.Bm,'LineWidth',lw,'Color',data.colorB); 
plot(t,s.Dm,'LineWidth',lw,'Color',data.colorD); 

legend('Reward','Higher-reward','Omission')
xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f);


%% Motor kernel
%{
s.rm2fit = s.rm(rng);
tfit = (0:length(s.rm2fit)-1)*dt;
[s.rkern, s.rgof] = f_alphafit(tfit,s.rm2fit);

figure
plot(t,s.rm,'Color',data.colorRandomlick,'LineWidth',lw); 
hold on;
plot(t,s.rkern(t-t0),'LineWidth',lw);
xlim(trange);

cont = 20;
figure
plot(t,l.DD(cont,:))
hold on
plot(t,s.r(cont,:),'Color','k')
%}

lick_fit = zeros(1,length(s.r(:,1)));

for i = 1:length(lick_fit)
   pos = find(t>0);
   aux = find(l.DD(i,pos)>0);
   if t(pos(1)+aux(1)-1) == 180 && t(pos(1)+aux(2)-1) > 360
       lick_fit(i) = 1;
   end
end

length(find(lick_fit >0))
lick_fit = logical(lick_fit);


figure
plot(t,mean(s.r(lick_fit,:),1))
hold on
stem(t,mean(l.DD(lick_fit,:),1))
xlim([-1000 1000]); 

s.rm_fit = mean(s.r(lick_fit,:),1);


s.rm2fit = s.rm_fit(rng);
tfit = (0:length(s.rm2fit)-1)*dt;
[s.rkern, s.rgof] = f_alphafit2(tfit,s.rm2fit);


t_aux = zeros(1,length(mean(l.DD(lick_fit,:),1)));
t_auxx = mean(l.DD(lick_fit,:),1)>0.8;

alpha = @(x)(max(0,s.rkern.a*(x-s.rkern.c)/s.rkern.b.*exp(1-(x-s.rkern.c)/s.rkern.b)));

alpha_conv = conv(alpha(t-t(1)+s.rkern.c),t_auxx);
length(alpha_conv)
length(t)
figure()
plot(alpha_conv)
alpha_conv = alpha_conv(1:length(t));
%plot(t+s.rkern.c+t0,alpha_conv)

figure
plot(t,s.rm_fit,'Color',data.colorRandomlick,'LineWidth',lw); 
hold on;
plot(t,s.rkern(t-t0),'LineWidth',lw);
xlim(trange);
plot(t,alpha(t-t0),'Color','g','LineWidth',lw)
hold on
plot(t,alpha(t-t0-180),'Color','g','LineWidth',lw)
stem(t,t_auxx)
plot(t+s.rkern.c+t0,alpha_conv,'Color','k','LineWidth',lw)

%% Convolution with motor kernel

s.motfitA = alpha(t-t(1)+s.rkern.c);
s.motfitB = alpha(t-t(1)+s.rkern.c);
s.motfitD = alpha(t-t(1)+s.rkern.c);

tc = t+s.rkern.c+t0;

s.motA = zeros(s.nA,s.npts);
s.motB = zeros(s.nB,s.npts);

for ii = 1:s.nA
    aux =  conv(l.A(ii,:),s.motfitA);
    s.motA(ii,:) = aux(1:length(l.A(1,:)));
end

figure
ax = s.motA(1,:);
plot(tc,ax)
hold on
aaa = l.A(1,:);
stem(t(aaa>0),aaa(aaa>0))

for ii = 1:s.nB
    aux =  conv(l.B(ii,:),s.motfitB);
    s.motB(ii,:) = aux(1:length(l.B(1,:)));
end

for ii = 1:s.nD
    aux =  conv(l.D(ii,:),s.motfitD);
    s.motD(ii,:) = aux(1:length(l.D(1,:)));
end

% To plot these right, we need to either shift time or the arrays. I shifted
% the arrays.
% tc = t+s.rkern.c+t0;
% index to shift by:
nshift = round((s.rkern.c+t0)/dt);
s.motA = circshift(s.motA,nshift,2); 
s.motB = circshift(s.motB,nshift,2);
s.motD = circshift(s.motD,nshift,2);
s.motAm = mean(s.motA,1);
s.motBm = mean(s.motB,1);
s.motDm = mean(s.motD,1);

figure;
plot(t,s.motAm,'LineWidth',lw,'Color',data.colorA);
hold on;
plot(t,s.motBm,'LineWidth',lw,'Color',data.colorB);
plot(t,s.motDm,'LineWidth',lw,'Color',data.colorD);
% Actual data:
plot(t,s.Am,'Color',data.colorA); 
plot(t,s.Bm,'Color',data.colorB); 
plot(t,s.Dm,'Color',data.colorD); 
hold off;
xlim(trange);  
legend('Reward','Higher-reward','Omission')
xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f);

%% Cue kernel    
s.cue2fit = s.Dm - s.motDm;

s.cue2fit = s.cue2fit(rng);
tfit = (0:length(s.cue2fit)-1)*dt;

[s.cue, s.cuegof] = f_alphafit(tfit,s.cue2fit);
s.cuefit = s.cue(t-t0)';

figure
plot(t,s.motAm+s.cuefit,'LineWidth',lw,'Color',data.colorA);
hold on;
plot(t,s.motBm+s.cuefit,'LineWidth',lw,'Color',data.colorB);
plot(t,s.motDm+s.cuefit,'LineWidth',lw,'Color',data.colorD);
% Actual data:
plot(t,s.Am,'Color',data.colorA); 
plot(t,s.Bm,'Color',data.colorB); 
plot(t,s.Dm,'Color',data.colorD); 
hold off;
xlim(trange);  
legend('Reward','Higher-reward','Omission')
xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f);

%% Reward signal

s.reward2fit = s.Am - s.motAm-s.cuefit;
s.hreward2fit = s.Bm - s.motBm-s.cuefit;

figure
plot(t,s.hreward2fit,'LineWidth',lw,'Color',data.colorB);
hold on
plot(t,s.reward2fit,'LineWidth',lw,'Color',data.colorA);

times = [];
n_licks = 2;
for i=1:length(l.A(:,1))
    aux = find(l.A(i,t>0)>0);
    t_aux = t(162+aux);
    times = [times, t_aux(1:n_licks)];
end

t0 = -150;
t1 = 4000;
rng = t>t0 & t<t1;
s.Am2fit = s.reward2fit(rng);
tfit = (0:length(s.Am2fit)-1)*dt;

[s.Akern, s.Agof] = f_alphafit3(tfit,s.Am2fit,times);

alpha_rew = @(x)(max(0,s.Akern.a*(x-s.Akern.c)/s.Akern.b.*exp(1-(x-s.Akern.c)/s.Akern.b)));


figure()
plot(t,s.motAm+s.cuefit+s.Akern(t-t0)','Color',data.colorA,'LineWidth',lw); 
hold on;
plot(t,s.motDm+s.cuefit,'LineWidth',lw,'Color',data.colorD);
plot(t,s.Am,'Color',data.colorA); 
plot(t,s.Bm,'Color',data.colorB); 
plot(t,s.Dm,'Color',data.colorD); 

xlim([-1000,3000]);  


%% Higher-Reward signal

s.hreward2fit = s.Bm - s.motBm-s.cuefit;

times = [];
for i=1:length(l.B(:,1))
    aux = find(l.B(i,t>0)>0);
    t_aux = t(162+aux);
    times = [times, t_aux(1:n_licks)];
end

t0h = -50;
t1h = t1;
rng = t>t0h & t<t1h;
s.Bm2fit = s.hreward2fit(rng);
tfit = (0:length(s.Bm2fit)-1)*dt;

[s.Bkern, s.Bgof] = f_alphafit3(tfit,s.Bm2fit,times);

alpha_hrew = @(x)(max(0,s.Bkern.a*(x-s.Bkern.c)/s.Bkern.b.*exp(1-(x-s.Bkern.c)/s.Bkern.b)));


figure()
plot(t,s.motAm+s.cuefit+s.Akern(t-t0)','Color',data.colorA,'LineWidth',lw); 
hold on
plot(t,s.motBm+s.cuefit+s.Bkern(t-t0h)','Color',data.colorB,'LineWidth',lw); 
plot(t,s.motDm+s.cuefit,'LineWidth',lw,'Color',data.colorD);
plot(t,s.Am,'Color',data.colorA); 
plot(t,s.Bm,'Color',data.colorB); 
plot(t,s.Dm,'Color',data.colorD); 

legend('Reward','Higher-reward','Omission','','','')
xlim(trange);  
ylabel('DCN-SNc');
xlabel('time[ms]')
set(gca,'fontsize',f);

figure
hold on
plot(t,alpha(t-t0),'Color','r','LineWidth',lw)
plot(t,s.cuefit,'Color','b','LineWidth',lw)
plot(t,10*alpha_rew(t),'Color',data.colorA,'LineWidth',lw);
plot(t,10*alpha_hrew(t),'Color',data.colorB,'LineWidth',lw);
legend('Motor','Cue','Reward','Higher-reward')
xlim(trange);  
xlabel('time[ms]')
set(gca,'fontsize',f);


