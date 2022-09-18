function y = rate_model_network_fit(x,tau2e,tau2i,tau3i,x2ethr,x2eslp,x2ithr,x2islp,x3ithr,x3islp)

J2i1e = 0;

J2e1e = 0.89;    % 0.9
J3i1e = 0.21;  % 0.5

J2i2e = 2.7;  %1
J2e2i = 4.3;  %1
J2e3i = 1.001;  % 0.63

Tmax = 10000;
dt = 0.1;

ton2ecue = 4000;
toff2ecue = ton2ecue + 250;
A2ecue = 0.04;

% Cue
t0_cue = 4000;
t1_cue = t0_cue+250;
J_cue  = 0.08;

% Reward / Higher reward
t0_rew = 4000+150;
t1_rew = t0_rew+250; 
J_rew = 0*0.08;%0*0.17;
J_hrew = 1*0.15;%0*0.26;

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
% Functions

heaviside=@(t) 0.5*(t == 0)+(t > 0);
pulse = @(t,t0,t1) heaviside(t-t0).*heaviside(t1-t);
relumax = @(x,m,k) min(1,max(0,0.5+(x-m)/(4*k)));
Halfsigexp=@(x,xth,xslp) 2*(1./(1+exp(-2*(x-xth)/xslp))-0.5).*(x-xth >=0);

t = 0:dt:Tmax;

% Inputs

%2e
F2ecue = A2ecue*(heaviside(t-ton2ecue).*heaviside(toff2ecue-t));

% Initial conditions
x2e = zeros(1,length(t));
x2i = zeros(1,length(t));
x3i = zeros(1,length(t));

x1e = zeros(1,length(t));
w = zeros(1,length(t));

x3i(1) = 0.3;
x2i(1) = 0.04;
x2e(1) = 0;

x1e(1) = 0;
w(1) = 0;


for j=1:length(t)-1
    
    k1x1e = (-x1e(j)+relumax(R11*x1e(j)+J_cue*pulse(t(j),t0_cue,t1_cue)+J_rew*pulse(t(j),t0_rew,t1_rew)+J_hrew*pulse(t(j),t0_rew,t1_rew)+R12*w(j),m1,k1) )/tau1;
    k1x2e = (-x2e(j)+Halfsigexp(J2e1e*x1e(j)-J2e2i*x2i(j)-J2e3i*x3i(j)+F2ecue(j),x2ethr,x2eslp))/tau2e;
    k1x2i = (-x2i(j)+Halfsigexp(J2i1e*x1e(j)+J2i2e*x2e(j),x2ithr,x2islp))/tau2i;
    k1x3i = (-x3i(j)+Halfsigexp(J3i1e*x1e(j),x3ithr,x3islp))/tau3i;
    k1w = (-w(j)+relumax(R21*x1e(j),m2,k2))/tau2;

    ax1e = x1e(j)+k1x1e*dt;
    ax2e = x2e(j)+k1x2e*dt;
    ax2i = x2i(j)+k1x2i*dt;
    ax3i = x3i(j)+k1x3i*dt;
    aw = w(j)+k1w*dt;

    k2x1e = (-ax1e+relumax(R11*ax1e+J_cue*pulse(t(j+1),t0_cue,t1_cue)+J_rew*pulse(t(j+1),t0_rew,t1_rew)+J_hrew*pulse(t(j+1),t0_rew,t1_rew)+R12*aw,m1,k1))/tau1; 
    k2x2e = (-ax2e+Halfsigexp(J2e1e*ax1e-J2e2i*ax2i-J2e3i*ax3i+F2ecue(j+1),x2ethr,x2eslp))/tau2e;
    k2x2i = (-ax2i+Halfsigexp(J2i1e*ax1e+J2i2e*ax2e,x2ithr,x2islp))/tau2i;
    k2x3i = (-ax3i+Halfsigexp(J3i1e*ax1e,x3ithr,x3islp))/tau3i;   
    k2w = (-aw+relumax(R21*ax1e,m2,k2))/tau2;   

    
    x1e(j+1) = x1e(j)+(k1x1e+k2x1e)*dt/2;
    x2e(j+1) = x2e(j)+(k1x2e+k2x2e)*dt/2;
    x2i(j+1) = x2i(j)+(k1x2i+k2x2i)*dt/2;
    x3i(j+1) = x3i(j)+(k1x3i+k2x3i)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;

end

% Output
aux2 = x2e(t>= 3000 &  t<=7000); 
y = aux2(1:600:(600*length(x)-1))';
