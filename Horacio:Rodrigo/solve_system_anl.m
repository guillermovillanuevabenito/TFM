function [Vanl,Vpeakanl,tpeakanl,Vlowanl,tlowanl ] = solve_system_anl(taurse, taudec, taudep, taufac, adep, afac, Xinf, Zinf, C, El, Gl, Iapp, Gsyn, Esyn, t,dt,tspk,Vpre,spkwidth,Tmax)

tspkpre = tspk;
% Analytical Calculation of M = DeltaS, X and Z       

ISIpre = diff(tspkpre(1:end-1));
X = zeros(1,length(tspkpre)-1);
Z = zeros(1,length(tspkpre)-1);
X(1) = 1;
Z(1) = afac;
for j=1:length(tspkpre)-2
    Dtatspkin = ISIpre(j);
    Pdep = exp(-Dtatspkin/taudep);
    Pfac =  exp(-Dtatspkin/taufac);
    X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*Pdep;       
    Z(j+1) = (1-afac)*(Zinf+(Z(j)-Zinf)*Pfac)+afac;
end
M = X.*Z;

%     plot(tspkpre(1:end-1),M,'o','Color',lightcoral,'linewidth',2)
%     legend('x','z','x z','x^- Z^+','M');
%     
% Analytical computation of S

%%%%%% Make M independent of the numerical calculation
   

jspkwidth = floor(spkwidth/dt);
%DeltaSspk = M*0.99;   
DeltaSspk = M*(taudec/(taudec+taurse)*(1-exp(-1*(taudec+taurse)/(taudec*taurse))));
Nspkeff = find(tspkpre>Tmax,1)-1;
jspk = floor(tspkpre/dt);
tt = 0:dt:Tmax;
Sanl = zeros(1,length(tt));
for j=1:jspk(1)
   Sanl(j) = 0;
end
for k=1:Nspkeff-1
    for j=jspk(k)+1:jspk(k)+jspkwidth
        Sanl(j) = DeltaSspk(k);
    end
    for j=jspk(k)+jspkwidth+1:jspk(k+1)
        Sanl(j) = DeltaSspk(k).*exp(-(tt(j)-tspkpre(k)-spkwidth)/taudec);
    end
end

% Analytical Computations of V, peaks and troughs  

eta = 1;

alpha = Gsyn/Gl*(Esyn-El);
tau = C/Gl;
Vanl = zeros(1,length(tt));
for j=1:jspk(1)
   Vanl(j) = 0;
end

Vo = zeros(1,length(tspkpre));
beta = zeros(1,length(tspkpre));
alphabase = alpha;
alpha = zeros(1); 
a = zeros(1);
b = zeros(1);
tpeakanl = zeros(1);
Vpeakanl = zeros(1);

Vlowanl = zeros(1);
tlowanl = zeros(1);

Vo(1) = 0; 
alpha(1) = alphabase;
sgma = 0;

cnt = 1;

if taudec ~= tau                
    for l=1:find(tspkpre<Tmax,1,'last')-1   
        alpha(l) = alphabase*(1-sgma);           
        for j=jspk(l)+1:jspk(l)+jspkwidth+1
            Vanl(j) = alpha(l)*DeltaSspk(l)+(Vo(l)-alpha(l)*DeltaSspk(l))*exp(-(tt(j)-tspkpre(l))/tau);
        end
        beta(l) = alpha(l)*DeltaSspk(l)+(Vo(l)-alpha(l)*DeltaSspk(l))*exp(-spkwidth/tau);
        for j=jspk(l)+jspkwidth+2:jspk(l+1)
            Vanl(j) = alpha(l)*taudec*DeltaSspk(l)/(taudec-tau)*exp(-(tt(j)-tspkpre(l)-spkwidth)/taudec)+(beta(l)-alpha(l)*taudec*DeltaSspk(l)/(taudec-tau))*exp(-(tt(j)-tspkpre(l)-spkwidth)/tau);
        end
        Vo(l+1) = lpaha(l)*taudec*DeltaSspk(l)/(taudec-tau)*exp(-(tspkpre(l+1)-tspkpre(l)-spkwidth)/taudec)+(beta(l)-alpha(l)*taudec*DeltaSspk(l)/(taudec-tau))*exp(-(tspkpre(l+1)-tspkpre(l)-spkwidth)/tau);                          
        a(l) = alpha(l)*taudec*DeltaSspk(l)/(taudec-tau)*exp(spkwidth/taudec);
        b(l) = (beta(l)-alpha(l)*taudec*DeltaSspk(l)/(taudec-tau))*exp(spkwidth/tau);
        % tpeakanl(l) = tspkpre(l)+taudec*tau/(taudec-tau)*log(-b(l)/tau*taudec/a(l));
        % Vpeakanl(l) = alpha(l)*taudec*DeltaSspk(l)/(taudec-tau)*exp(-(tpeakanl(l)-tspkpre(l)-spkwidth)/taudec)+(beta(l)-alpha(l)*taudec*DeltaSspk(l)/(taudec-tau))*exp(-(tpeakanl(l)-tspkpre(l)-spkwidth)/tau);                                
        %sgma = ((1-eta)*Vo(l+1)+eta*Vpeakanl(l))/(Esyn-El);  
        [Vpeakanl(l),jpeak] = max(Vanl(jspk(l)+1:jspk(l+1)));
        Vlowanl(l) = Vanl(jspk(l+1));
        tpeakanl(l) = t(jspk(l)+1+jpeak);
        tlowanl(l) = t(jspk(l+1));
        sgma = (eta*Vpeakanl(l))/(Esyn-El);     
        Vtroughanl(l) = Vo(l+1); 
    end
 else        
    for l=1:find(tspkpre<Tmax,1,'last')-1   
        alpha(l) = alphabase*(1-sgma);            
        for j=jspk(l)+1:jspk(l)+jspkwidth+1
            Vanl(j) = alpha(l)*DeltaSspk(l)+(Vo(l)-alpha(l)*DeltaSspk(l))*exp(-(tt(j)-tspkpre(l))/tau);
        end
        beta(l) = alpha(l)*DeltaSspk(l)+(Vo(l)-alpha(l)*DeltaSspk(l))*exp(-spkwidth/tau);
        for j=jspk(l)+jspkwidth+2:jspk(l+1)                
            Vanl(j) = alpha(l)*DeltaSspk(l)/tau*exp(spkwidth/taudec)*tt(j)*exp(-(tt(j)-tspkpre(l))/tau)+(beta(l)-alpha(l)*DeltaSspk(l)/tau*(tspkpre(l)+spkwidth))*exp(spkwidth/tau)*exp(-(tt(j)-tspkpre(l))/tau);                           
        end
        Vo(l+1) = alpha(l)*DeltaSspk(l)/tau*exp(spkwidth/taudec)*tspkpre(l+1)*exp(-(tspkpre(l+1)-tspkpre(l))/tau)+(beta(l)-alpha(l)*DeltaSspk(l)/tau*(tspkpre(l)+spkwidth))*exp(spkwidth/tau)*exp(-(tspkpre(l+1)-tspkpre(l))/tau);                          
        a(l) = alpha(l)*DeltaSspk(l)/tau*exp(spkwidth/taudec);
        b(l) = (beta(l)-alpha(l)*DeltaSspk(l)/tau*(tspkpre(l)+spkwidth))*exp(spkwidth/tau);
        % tpeakanl(l) = tau-b(l)/a(l);
        % Vpeakanl(l) = alpha(l)*DeltaSspk(l)/tau*exp(spkwidth/taudec)*tpeakanl(l)*exp(-(tpeakanl(l)-tspkpre(l))/tau)+(beta(l)-alpha(l)*DeltaSspk(l)/tau*(tspkpre(l)+spkwidth))*exp(spkwidth/tau)*exp(-(tpeakanl(l)-tspkpre(l))/tau);                                                                            
        %sgma = ((1-eta)*Vo(l+1)+eta*Vpeakanl(l))/(Esyn-El);     


        [aux1,jpeak] = max(Vanl(jspk(l)+1:jspk(l+1)));
        aux2 = Vanl(jspk(l+1));

        if aux1 ~= aux2
            Vpeakanl(cnt) = aux1;
            Vlowanl(cnt) = aux2;
            tpeakanl(cnt) = t(jspk(l)+1+jpeak);
            tlowanl(cnt) = t(jspk(l+1));


            sgma = (eta*Vpeakanl(cnt))/(Esyn-El);             
            Vtroughanl(cnt) = Vo(l+1); 

            cnt = cnt+1;
        end
     

        %[Vpeakanl(l),jpeak] = max(Vanl(jspk(l)+1:jspk(l+1)));
        %Vlowanl(l) = Vanl(jspk(l+1));
        %tpeakanl(l) = t(jspk(l)+1+jpeak);
        %tlowanl(l) = t(jspk(l+1));

        %sgma = (eta*Vpeakanl(cnt))/(Esyn-El);             
        %Vtroughanl(cnt) = Vo(l+1);
        a(l) = alpha(l)*DeltaSspk(l)/tau*exp((tspkpre(l)+spkwidth)/taudec);
        b(l) = (beta(l)-alpha(l)*DeltaSspk(l)/tau*(tspkpre(l)+spkwidth))*exp((tspkpre(l)+spkwidth)/tau);
    end

end
