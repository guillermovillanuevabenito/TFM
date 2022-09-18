% solve the model
function [V] = solve_system(taurse, taudec, taudep, taufac, adep, afac, Xinf, Zinf, C, El, Gl, Iapp, Gsyn, Esyn, t,dt,tspk,Vpre)

S = zeros(1,length(t));
x = zeros(1,length(t));
z = zeros(1,length(t));
M = zeros(1,length(tspk));
V = zeros(1,length(t));


S(1) = 0;
x(1) = 1;
z(1) = 0;
V(1) = El;

DeltaS = 0;
jspk = 1;


% Functions

    % Sigmoid function

vsth = 0;
vsslp = 4;
H=@(v) 0.5*(1 + tanh(v-vsth/vsslp));

for j=1:length(t)-1
    k1x = (Xinf-x(j))/taudep;
    k1z = (Zinf-z(j))/taufac;
    k1s = H(Vpre(j))*(DeltaS-S(j))/taurse-S(j)/taudec;
    k1v = (-Gl*(V(j)-El)+Iapp-Gsyn*S(j)*(V(j)-Esyn))/C;
    as = S(j)+k1s*dt;
    ax = x(j)+k1x*dt;
    az = z(j)+k1z*dt;
    av = V(j)+k1v*dt;
    k2x = (Xinf-ax)/taudep;
    k2z = (Zinf-az)/taufac;
    k2s = H(Vpre(j))*(DeltaS-as)/taurse-as/taudec; 
    k2v = (-Gl*(av-El)+Iapp-Gsyn*as*(av-Esyn))/C;
    S(j+1) = S(j)+(k1s+k2s)*dt/2;  
    x(j+1) = x(j)+(k1x+k2x)*dt/2;
    z(j+1) = z(j)+(k1z+k2z)*dt/2;
    V(j+1) = V(j)+(k1v+k2v)*dt/2;  
    if jspk < length(tspk) && t(j)>= tspk(jspk)
        x(j+1) = x(j+1)-adep*x(j);
        z(j+1) = z(j+1)+afac*(1-z(j+1));            
        DeltaS = x(j)*z(j+1);
        M(jspk) = DeltaS;
        jspk = jspk+1;
    end
end