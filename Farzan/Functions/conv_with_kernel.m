function y = conv_with_kernel(x,a,b,c)
aux = [0,180];
y = 0*x;
for i=1:length(aux)
     y = y + max(0,a*(x-c-aux(i))/b.*exp(1-(x-c-aux(i))/b) );
end
%y = aux1;
%y = max(0,a*(x-c-aux(1))/b.*exp(1-(x-c-aux(1))/b) ) + max(0,a*(x-c-aux(2))/b.*exp(1-(x-c-aux(2))/b) )
