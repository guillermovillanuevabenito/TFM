% Given a voltage trace, computes peak and low times and their values
function [Vpeak, tpeak, Vlow, tlow] = compute_metrics(V,t)

El = -60;

tpeak = zeros(1);
Vpeak = zeros(1);

tlow = zeros(1);
Vlow = zeros(1);

cntpeak = 0;

for j=1:length(t)-1
        
    if t(j)>2  

        if V(j)>V(j-1) && V(j)>V(j+1)

            cntpeak = cntpeak+1;
            tpeak(cntpeak) = t(j);
            Vpeak(cntpeak) = V(j);  
        end 
        if cntpeak > 0
            if (V(j)< V(j-1) && V(j) < V(j+1)) || abs(V(j)-El) < 1e-3
                tlow(cntpeak) = t(j);
                Vlow(cntpeak) = V(j);  
            end
        end

    end
end

if length(Vpeak) ~= length(Vlow)
    Vpeak = Vpeak(1:end-1);
    cntpeak = cntpeak-1;
    tpeak = tpeak(1:end-1);
end