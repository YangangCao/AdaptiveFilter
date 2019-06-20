%
% SQUARE-ROOT RAISED-COSINE PULSE: h=sr_cos_p(N,L,alpha)
% This function generates a square-root raised-cosine pulse of length N+1. 
% There are L samples per symbol period.
% alpha is the roll-off factor.
%
function h=sr_cos_p(N,L,alpha)
t=[-N/2:1:N/2]/L;
h=zeros(size(t));
for k=1:length(t)
    if t(k)==0
        h(k)=1-alpha+4*alpha/pi;
    elseif (t(k)==(-1/(4*alpha)))|(t(k)==(1/(4*alpha)))
        h(k)=(alpha/sqrt(2))*((1+2/pi)*sin(pi/4/alpha)+(1-2/pi)*cos(pi/4/alpha));
    else
        h(k)=(sin(pi*(1-alpha)*t(k))+4*alpha*t(k)*cos(pi*(1+alpha)*t(k)))/(pi*t(k)*(1-(4*alpha*t(k))^2));
    end
end
h=h'/sqrt(L);