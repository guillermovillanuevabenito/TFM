function y = rate_fit(x,tau1,tau2,m1,k1,m2,k2)

%%%%% Parameters %%%%%

t_start = 4000;
t_sim = 4000;

t_max = t_start+t_sim;
dt = 0.1;

% Cue
t0_cue = t_start;
t1_cue = t0_cue+250;
J_cue  = 0.2;            % good 0.2;

% Reward / Higher reward
t0_rew = t_start+150;
t1_rew = t0_rew+250; 
J_rew = 0*0.1;
J_hrew = 0*0.12;                % good 0*0.23;

% Activation function DCN
%m1 = 0.34;
%k1 = 0.15;
%tau = 0;

% Activation function recovery variable
%m2 = 0.215;
%k2 = 0.12;
%tau2 = 500;

% Connectivities
R11 = 4;
R12 = -1.8;
R21 = 1;

t= 0:dt:t_max;

% Functions
pulse = @(t,t0,t1) heaviside(t-t0).*heaviside(t1-t);
relumax = @(x,m,k) min(1,max(0,0.5+(x-m)/(4*k)));

% Initial conditions
x1e = zeros(1,length(t));
w = zeros(1,length(t));

x1e(1) = 0;
w(1) = 0;

for j=1:length(t)-1
   
    k1x1e = (-x1e(j)+relumax(R11*x1e(j)+J_cue*pulse(t(j),t0_cue,t1_cue)+J_rew*pulse(t(j),t0_rew,t1_rew)+J_hrew*pulse(t(j),t0_rew,t1_rew)+R12*w(j),m1,k1) )/tau1;
    k1w = (-w(j)+relumax(R21*x1e(j),m2,k2))/tau2;
    
    ax1e = x1e(j)+k1x1e*dt;
    aw = w(j)+k1w*dt;
    
    k2x1e = (-ax1e+relumax(R11*ax1e+J_cue*pulse(t(j+1),t0_cue,t1_cue)+J_rew*pulse(t(j+1),t0_rew,t1_rew)+J_hrew*pulse(t(j+1),t0_rew,t1_rew)+R12*aw,m1,k1))/tau1; 
    k2w = (-aw+relumax(R21*ax1e,m2,k2))/tau2;   
    
    x1e(j+1) = x1e(j)+(k1x1e+k2x1e)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
        
end

% Output
aux = x1e(t>= 3000 &  t<=7000); 
y = aux(1:600:(600*length(x)-1))';
%y = aux(1:floor(length(aux)/60):(floor(length(aux)/60)*(length(x)-1)+1) )';
