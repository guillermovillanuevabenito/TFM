function [y,y1,y2] = rate_model_network(x,DCN,Tmax,dt,x2ethr,x2eslp,x2ithr,x2islp,x3ithr,x3islp,tau2e,tau2i,tau3i,J2i1e,J2e1e,J2i2e,J2e2i,J2e3i,J3i1e,ton2ecue,toff2ecue,A2ecue)

% Functions

heaviside=@(t) 0.5*(t == 0)+(t > 0);
Halfsigexp=@(x,xth,xslp) 2*(1./(1+exp(-2*(x-xth)/xslp))-0.5).*(x-xth >=0);
relumax = @(x,m,k) min(1,max(0,0.5+(x-m)/(4*k)));


t = 0:dt:Tmax;

% Inputs

%2e
F2ecue = A2ecue*(heaviside(t-ton2ecue).*heaviside(toff2ecue-t));

% Initial conditions
x2e = zeros(1,length(t));
x2i = zeros(1,length(t));
x3i = zeros(1,length(t));

x3i(1) = 0.3;
x2i(1) = 0.04;
x2e(1) = 0;


for j=1:length(t)-1
    
    k1x2e = (-x2e(j)+Halfsigexp(J2e1e*DCN(j)-J2e2i*x2i(j)-J2e3i*x3i(j)+F2ecue(j),x2ethr,x2eslp))/tau2e;
    k1x2i = (-x2i(j)+Halfsigexp(J2i1e*DCN(j)+J2i2e*x2e(j),x2ithr,x2islp))/tau2i;
    k1x3i = (-x3i(j)+Halfsigexp(J3i1e*DCN(j),x3ithr,x3islp))/tau3i;

    ax2e = x2e(j)+k1x2e*dt;
    ax2i = x2i(j)+k1x2i*dt;
    ax3i = x3i(j)+k1x3i*dt;
    
    k2x2e = (-ax2e+Halfsigexp(J2e1e*DCN(j+1)-J2e2i*ax2i-J2e3i*ax3i+F2ecue(j+1),x2ethr,x2eslp))/tau2e;
    k2x2i = (-ax2i+Halfsigexp(J2i1e*DCN(j+1)+J2i2e*ax2e,x2ithr,x2islp))/tau2i;
    k2x3i = (-ax3i+Halfsigexp(J3i1e*DCN(j+1),x3ithr,x3islp))/tau3i;   

        
    x2e(j+1) = x2e(j)+(k1x2e+k2x2e)*dt/2;
    x2i(j+1) = x2i(j)+(k1x2i+k2x2i)*dt/2;
    x3i(j+1) = x3i(j)+(k1x3i+k2x3i)*dt/2;

end

% Output
aux = x2e(t>= 3000 &  t<=7000); 
y = aux(1:600:(600*length(x)-1))';

% Output
aux1 = x2i(t>= 3000 &  t<=7000); 
y1 = aux1(1:600:(600*length(x)-1))';

% Output
aux2 = x3i(t>= 3000 &  t<=7000); 
y2 = aux2(1:600:(600*length(x)-1))';