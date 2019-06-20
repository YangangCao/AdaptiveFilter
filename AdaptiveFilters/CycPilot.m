%
% This function designs a cyclic pilot sequence of length N+1.
% It follows the construction formula (11.55).
% 
function s=CycPilot(N)
if rem(N+1,2)==0
    s=exp(j*pi*([0:N]').^2/(N+1));
else
    s=exp(j*pi*[0:N]'.*[1:N+1]'/(N+1));
end
