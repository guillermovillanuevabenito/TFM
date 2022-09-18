function [signal, d_signal, rec, t] = rate_model_2d(t_max, dt, t0_cue, t1_cue, J_cue,t0_rew, t1_rew, J_rew, J_hrew, m1, k1, tau1, m2,k2,tau2, R11, R12, R21)

t = 0:dt:t_max;

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

D_x1e = (-x1e+relumax(R11*x1e+J_cue*pulse(t,t0_cue,t1_cue)+J_rew*pulse(t,t0_rew,t1_rew)+J_hrew*pulse(t,t0_rew,t1_rew)+R12*w,m1,k1) )/tau1;

% Output
signal = x1e;
d_signal = D_x1e;
rec = w;
