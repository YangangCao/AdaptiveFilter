function h=rcos50(t);
%
% This function gives samples of a raised cosine 
% pulse shape with a rolloff factor of 50%.
% 
% The sampling times are specified by the elements
% of the vector "t".
%
t=t+eps;
h=(sin(pi*t)./(pi*t)).*(cos(pi*t/2)./(1-t.^2));

