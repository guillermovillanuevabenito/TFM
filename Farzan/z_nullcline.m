function [z] = z_nullcline(xx,m1,k1,m2,k2,R11,R12,R21,I)
z = 0.*xx;

for i=1:length(xx)
    z(i) = xx(i)*R12/(4*k1-R11) + (2*k1-m1+I)/(4*k1-R11);
end


for i=1:length(z)
    if z(i)> (2*k1+m1-R11)/R12
        z(i) = 1;
    end
end