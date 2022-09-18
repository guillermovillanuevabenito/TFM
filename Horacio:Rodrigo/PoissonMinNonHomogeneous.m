% Compute poisson spike inputs
function [tspk,ISI] = PoissonMinNonHomogeneous(rconst,rmod,Tmax,dt,ISImin,R)

t = 0:dt:Tmax;
rconst = rconst/(length(t)-1)*Tmax;         % [rate] = Spk/bin
rmod = rmod/(length(t)-1)*Tmax;             % [rate] = Spk/bin
r = rconst+rmod*R;
spk = double(rand(1,length(t))<r);          % Generated spike trains
tspkbase = find(spk'>0)*dt;                 % Spike times
ISIbase = diff(tspkbase);
ISI = zeros(1);
l=0;
for j=1:length(ISIbase)
    if ISIbase(j) > ISImin
        l=l+1;
        ISI(l) = ISIbase(j);
    end
end    
tspk = zeros(1);
tspk(1) = tspkbase(1);
for l=2:length(ISI)
    tspk(l) = tspk(l-1)+ISI(l);
end



