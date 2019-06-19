function [W, e] = lms(u, d, mu);
% Maximum number of time step that can be predicted
N = min(size(u, 1),size(d, 1));
Nin = size(u,2);
Nout =size(d,2);
% Intializatize weight matrix and associated parameters for LMS predictor
w = zeros(Nout, Nin);
W = [];
for n=1:N,
W = [W;w];
% Predict next sample and error
xp(n, :) = u(n,:)*w';
e(n,:) = d(n,:)-xp(n,:);
% Adapt weight matrix ans step size
w = w + mu * e(n,:)' * u(n,:);
end; % for n t
return; 
