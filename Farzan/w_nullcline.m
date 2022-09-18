function [z] = w_nullcline(xx,m1,k1,m2,k2,R11,R12,R21)
z = 0.*xx;
for i=1:length(xx)
    if xx(i) < (m2-2*k2)
        z(i) = 0;
    elseif xx(i) > (2*k2+m2)
        z(i) = 1;
    else
        z(i) = 0.5 + (R21*xx(i)-m2)/(4*k2);
    end
end