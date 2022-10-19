
a = 1;
b = 250;
c = -50;

alpha_function = @(x)(max(0,a*(x-c)/b.*exp(1-(x-c)/b)));

t = -1000:0.01:3000;
figure(1)
hold on
plot(t,10*alpha_function(t),'k','LineWidth',lw);
xlabel('time[ms]')
set(gca,'fontsize',f);
