clearvars;
close all;


% Functions

    % Heaviside
    
heaviside=@(t) 0.5*(t == 0)+(t > 0);

    % Sigmoid function

vsth = 0;
vsslp = 4;
H=@(v) 0.5*(1 + tanh(v-vsth/vsslp));

    % Depression and facilitation (DA model)
    
Q=@(a,tau,Delta) (1-a).*exp(-Delta/tau);
Xbar=@(a,tau,Delta,Xinf) (1-exp(-Delta/tau))*Xinf./(1-(1-a).*exp(-Delta/tau));
Zbar=@(a,tau,Delta,Zinf) ((1-exp(-Delta/tau))*(1-a)*Zinf+a)./(1-(1-a).*exp(-Delta/tau));

taudep_range = [100];

max_amp_num_var = [];
max_peak_num_var = [];
max_amp_num_mean = [];

max_amp_jit_var = [];
max_amp_jit_mean = [];

sigma_range = 0.25;


var_mean_jitter_amp = [];
var_std_jitter_amp = [];

var_mean_jitter_peak = [];
var_std_jitter_peak = [];

%%%
for sigma=sigma_range

sigma

%taufac = taudep;
v_peak_anl_mean = [];
v_peak_anl_var = [];

v_peak_num_mean = [];
v_peak_num_var = [];

v_amp_num_mean = [];
v_amp_num_var = [];

v_amp_anl_mean = [];
v_amp_anl_var = [];

v_peak_per_mean = [];
v_amp_per_mean = [];

v_peak_jit_mean = [];
v_peak_jit_var = [];

v_amp_jit_mean = [];
v_amp_jit_var = [];

%%%

SpkFreqin_range = 1:10:200;
%SpkFreqin_range = 100;


%%% Loop in frequencies
for SpkFreqin = SpkFreqin_range

    %SpkFreqin

    % Parameters
     
    taurse = 0.1;
    taudec = 10;

    taufac = 100;
    taudep = taufac;

    adep = 0.1;
    afac = 0.1;
    Xinf = 1;
    Zinf = 0; 

    C = 1;
    El = -60;
    Gl = 0.1;
    Iapp = 0;
    Gsyn = 0.1;
    Esyn = 0;
   
    % Time definitions

    Tmax = 100000;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % Input spike train

    spkwidth = 1;
    Tolerance = 0.001; 
    
    % Generation of Poisson-distributed spike trains        
    
    rconst = SpkFreqin/1000;      
    rmod = 0;
    R = zeros(1,length(t));
  
    ISImin = 5;
    [tspkpre,ISIpre] = PoissonMinNonHomogeneous(rconst,rmod,Tmax,dt,ISImin,R); 
    Vspkpre = zeros(1,length(t));
    spkprewidth = spkwidth;    

    for k=1:length(tspkpre)
        Vspkpre(floor(tspkpre(k)/dt):floor(tspkpre(k)/dt)+floor(spkprewidth/dt))=1;
    end
    Vspkpre = 120*Vspkpre-60;

    % Numerical computations: periodic spike inputs
    
    Tmax = 98000;    
    dt = 0.01;
    t = 0:dt:Tmax;
    
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspk = (1:Nspk+1)*SpkPerin-SpkPerin/2; 
    Vpre = 60*square(2*pi*SpkFreqin*t/1000-pi,100*spkwidth/SpkPerin);  

    %{
    figure
    hold on
    plot(t,Vpre,'-b','linewidth',2);
    plot(tspk,60*ones(length(tspk)),'or','linewidth',2)

    axis([0 100 0 60])
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    %}
    
    V = solve_system(taurse, taudec, taudep, taufac, adep, afac, Xinf, Zinf, C, El, Gl, Iapp, Gsyn, Esyn, t,dt,tspk, Vpre);
    [Vpeak_per, tpeak_per, Vlow_per, tlow_per] = compute_metrics(V,t);

    %{
    figure
    hold on
    plot(t,V,'-b','linewidth',2);
    hold on
    plot(tpeak_per,Vpeak_per,'ro','linewidth',2)
    plot(tlow_per,Vlow_per,'go','linewidth',2)

    axis([0 Tmax -70 -40])
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    %}


    % Numerical computations: jittered spike inputs
    tspk_jit = tspk; 
    tspk_jit(2:end) = tspk_jit(2:end)+randn(1,length(tspk_jit(2:end)))*(sigma)*SpkPerin;
    Vpre_jit = zeros(1,length(t));

    for k=1:length(tspk)
        Vpre_jit(floor(tspk_jit(k)/dt):floor(tspk_jit(k)/dt)+floor(spkwidth/dt))=1;
    end
    Vpre_jit = Vpre_jit*120-60;

    %{
    figure
    hold on
    plot(t,Vpre_jit(1:length(t)),'-b','linewidth',2);
    plot(tspk_jit,60*ones(length(tspk_jit)),'or','linewidth',2)


    axis([0 100 0 70])
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    %}

    V_jit = solve_system(taurse, taudec, taudep, taufac, adep, afac, Xinf, Zinf, C, El, Gl, Iapp, Gsyn, Esyn, t,dt,tspk_jit, Vpre_jit);
    [Vpeak_jit, tpeak_jit, Vlow_jit, tlow_jit] = compute_metrics(V_jit,t);

    %{
    figure
    hold on
    plot(t,V_jit,'-b','linewidth',2);
    hold on
    plot(tpeak_jit,Vpeak_jit,'ro','linewidth',2)
    plot(tlow_jit,Vlow_jit,'go','linewidth',2)

    axis([0 Tmax -70 -40])
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    %}

    %%% Compute poisson spike inputs, numerical solution
    V_pos_num = solve_system(taurse, taudec, taudep, taufac, adep, afac, Xinf, Zinf, C, El, Gl, Iapp, Gsyn, Esyn, t,dt,tspkpre, Vspkpre);
    [Vpeak_pos_num, tpeak_pos_num, Vlow_pos_num, tlow_pos_num] = compute_metrics(V_pos_num,t);


    %%% Compute poisson spike inputs, analytical solution
    %[V_pos_anl,Vpeak_pos_anl, tpeak_pos_anl, Vlow_pos_anl, tlow_pos_anl] = solve_system_anl(taurse, taudec, taudep, taufac, adep, afac, Xinf, Zinf, C, El, Gl, Iapp, Gsyn, Esyn, t,dt,tspkpre, Vspkpre,spkwidth,Tmax);
    %[Vpeak_pos_anl2, tpeak_pos_anl2, Vlow_pos_anl2, tlow_pos_anl2] = compute_metrics(V_pos_anl,t);
    
   %{
    figure
    hold on
    plot(t,V_pos_num,'-b','linewidth',2);
    hold on
    plot(tpeak_pos_num,Vpeak_pos_num,'ro','linewidth',2,'MarkerSize',8)
    plot(tlow_pos_num,Vlow_pos_num,'go','linewidth',2,'MarkerSize',8)
    %plot(tpeak_pos_anl2,Vpeak_pos_anl2-60,'k*','linewidth',2)
    %plot(tlow_pos_anl2,Vlow_pos_anl2-60,'m*','linewidth',2)


    axis([3000 3300 -57 -45])
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');

    figure
    hold on
    plot(t,V,'-b','linewidth',2);
    hold on
    plot(tpeak_per,Vpeak_per,'ro','linewidth',2,'MarkerSize',8)
    plot(tlow_per,Vlow_per,'go','linewidth',2,'MarkerSize',8)
    %plot(tpeak_pos_anl2,Vpeak_pos_anl2-60,'k*','linewidth',2)
    %plot(tlow_pos_anl2,Vlow_pos_anl2-60,'m*','linewidth',2)


    axis([3000 3300 -53 -48])
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    %}
    
%}
    tend = 3000;

    v_peak_per_mean = [v_peak_per_mean, mean(Vpeak_per(tpeak_per > tend))];
    aux = Vpeak_per-Vlow_per;
    v_amp_per_mean = [v_amp_per_mean, mean(aux(tpeak_per>tend)) ];
    
    v_peak_jit_mean = [v_peak_jit_mean, mean(Vpeak_jit(tpeak_jit > tend))];
    v_peak_jit_var = [v_peak_jit_var, var(Vpeak_jit(tpeak_jit > tend))];
    
    aux = Vpeak_jit - Vlow_jit;
    v_amp_jit_mean = [v_amp_jit_mean, mean(aux(tpeak_jit>tend))];
    v_amp_jit_var = [v_amp_jit_var,var(aux(tpeak_jit>tend))];


    v_peak_num_mean = [v_peak_num_mean, mean(Vpeak_pos_num(tpeak_pos_num > tend)) ];
    v_peak_num_var = [v_peak_num_var, var(Vpeak_pos_num(tpeak_pos_num > tend))];
    
    %v_peak_anl_mean = [v_peak_anl_mean, mean(Vpeak_pos_anl(tpeak_pos_anl > tend)+El) ];
    %v_peak_anl_var = [v_peak_anl_var, var(Vpeak_pos_anl(tpeak_pos_anl > tend)+El)];

    aux = Vpeak_pos_num-Vlow_pos_num;
    v_amp_num_mean = [v_amp_num_mean, mean(aux(tpeak_pos_num>tend))];
    v_amp_num_var = [v_amp_num_var, var(aux(tpeak_pos_num>tend))];
    
    %aux = Vpeak_pos_anl-Vlow_pos_anl;
    %v_amp_anl_mean = [v_amp_anl_mean, mean(aux(tpeak_pos_anl>tend))];
    %v_amp_anl_var = [v_amp_anl_var, var(aux(tpeak_pos_anl>tend))];


end


f1 = figure();
f1.Position = [10 10 1000 400];
subplot(1,2,1);
x = SpkFreqin_range;
plot(x,v_amp_per_mean,'g', 'LineWidth', 2)
hold on
plot(x, v_amp_jit_mean, 'k', 'LineWidth', 2);
plot(x, v_amp_num_mean, 'b', 'LineWidth', 2);
%plot(x, v_amp_anl_mean, 'r', 'LineWidth', 2);


ylabel('V,amp');
xlabel('Rate [Hz]')
set(gca,'fontsize',24);

%curve1 = v_amp_anl_mean-sqrt(v_amp_anl_var);
%curve2 = v_amp_anl_mean+sqrt(v_amp_anl_var);
%x2 = [x, fliplr(x)];
%inBetween = [curve1, fliplr(curve2)];
%fill(x2, inBetween, 'r','FaceAlpha',0.1,'LineStyle','none');


curve1 = v_amp_num_mean-sqrt(v_amp_num_var);
curve2 = v_amp_num_mean+sqrt(v_amp_num_var);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b','FaceAlpha',0.1,'LineStyle','none');


curve1 = v_amp_jit_mean-sqrt(v_amp_jit_var);
curve2 = v_amp_jit_mean+sqrt(v_amp_jit_var);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.2,'LineStyle','none');
legend('Periodic','Jittered','Poisson','','')

xlabel('Rate [Hz]')
set(gca,'fontsize',24);



subplot(1,2,2);
x = SpkFreqin_range;
plot(x, v_peak_per_mean, 'g', 'LineWidth', 2);
hold on
plot(x,v_peak_jit_mean,'k','LineWidth', 2)
plot(x,v_peak_num_mean,'b','LineWidth', 2)
%plot(x,v_peak_anl_mean,'r','LineWidth', 2)


ylabel('V,peak');
xlabel('Rate [Hz]')
set(gca,'fontsize',24);

%curve1 = v_peak_anl_mean-sqrt(v_peak_anl_var);
%curve2 = v_peak_anl_mean+sqrt(v_peak_anl_var);
%x2 = [x, fliplr(x)];
%inBetween = [curve1, fliplr(curve2)];
%fill(x2, inBetween, 'r','FaceAlpha',0.1,'LineStyle','none');


curve1 = v_peak_num_mean-sqrt(v_peak_num_var);
curve2 = v_peak_num_mean+sqrt(v_peak_num_var);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b','FaceAlpha',0.1,'LineStyle','none');


curve1 = v_peak_jit_mean-sqrt(v_peak_jit_var);
curve2 = v_peak_jit_mean+sqrt(v_peak_jit_var);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.2,'LineStyle','none');

legend('Periodic','Jittered','Poisson','','')
xlabel('Rate [Hz]')
set(gca,'fontsize',24);
%sgtitle(['taudep = ',num2str(taudep),'; taufac = ',num2str(taufac), '; Jittered: (std = 0.25*SpkPerin)'],'fontsize',24)


max_amp_num_var = [max_amp_num_var,SpkFreqin_range(find(v_amp_num_var==max(v_amp_num_var)))];
max_peak_num_var = [max_peak_num_var,SpkFreqin_range(find(v_peak_num_var==max(v_peak_num_var)))];
max_amp_num_mean = [max_amp_num_mean,SpkFreqin_range(find(v_amp_num_mean==max(v_amp_num_mean)))];

max_amp_jit_var = [max_amp_jit_var,SpkFreqin_range(find(v_amp_jit_var==max(v_amp_jit_var)))];
max_amp_jit_mean = [max_amp_jit_mean,SpkFreqin_range(find(v_amp_jit_mean==max(v_amp_jit_mean)))];


var_mean_jitter_amp = [var_mean_jitter_amp, mean(v_amp_jit_var)];
var_std_jitter_amp = [var_std_jitter_amp, sqrt(var(v_amp_jit_var))];

var_mean_jitter_peak = [var_mean_jitter_peak, mean(v_peak_jit_var)];
var_std_jitter_peak = [var_std_jitter_peak, sqrt(var(v_peak_jit_var))];

end

%{
figure
plot(taudep_range, max_amp_num_var, '-ob', 'LineWidth', 2);
hold on
plot(taudep_range, max_amp_jit_var, '-ok', 'LineWidth', 2);

legend('Poisson','Jittered')
ylabel('Rate [Hz]');
xlabel('\tau_{dep}=\tau_{fac}')
set(gca,'fontsize',20)

figure
plot(taudep_range, max_peak_num_var, 'g', 'LineWidth', 2);
hold on

figure
plot(taudep_range, max_amp_num_mean, '-ob', 'LineWidth', 2);
hold on
plot(taudep_range, max_amp_jit_mean, '-ok', 'LineWidth', 2);
legend('Poisson','Jittered')
ylabel('Rate [Hz]');
xlabel('\tau_{dep}=\tau_{fac}')
set(gca,'fontsize',20)

%}

%{

figure
plot(sigma_range*SpkPerin, var_mean_jitter_amp, '-ok', 'LineWidth', 2);
hold on
xlabel('\sigma')
set(gca,'fontsize',20)

figure
hold on
plot(sigma_range*SpkPerin, var_mean_jitter_peak, '-ok', 'LineWidth', 2);
xlabel('\sigma/\Delta_{spk}')
set(gca,'fontsize',20)



f1 = figure();
f1.Position = [10 10 1000 400];
subplot(1,2,1);
x = sigma_range;
plot(x,var_mean_jitter_amp,'k', 'LineWidth', 2)
hold on

xlabel('\sigma/\Delta_{spk}')
ylabel('Amp Var')

set(gca,'fontsize',24);

curve1 = var_mean_jitter_amp-var_std_jitter_amp;
curve2 = var_mean_jitter_amp+var_std_jitter_amp;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.1,'LineStyle','none');

set(gca,'fontsize',24);



subplot(1,2,2);
plot(x, var_mean_jitter_peak, 'k', 'LineWidth', 2);
hold on
xlabel('\sigma/\Delta_{spk}')
ylabel('Peak Var')

set(gca,'fontsize',24);

curve1 = var_mean_jitter_peak-var_std_jitter_peak;
curve2 = var_mean_jitter_peak+var_std_jitter_peak;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.1,'LineStyle','none');

set(gca,'fontsize',24);

%}