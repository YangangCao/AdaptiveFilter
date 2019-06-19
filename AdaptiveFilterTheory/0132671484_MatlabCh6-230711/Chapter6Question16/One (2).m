% lms_AR_pred.m - use multidimensional LMS algorithm to predict AR process
% written for MATLAB 4.0
%
% Input parameters:
% Xi : matrix of training/test points - each row is
% considered a datum
% Xi0 : initial value Xi(0) of series
% verbose : set to 1 for interactive processing
% mu : the initial value of step size
% decay : set to 1 for O(1/n) decay in mu
Nout = size(Xi, 2);
% length of maximum number of timesteps that can be predicted
N = size(Xi, 1);
% initialize weight matrix and associated parameters for LMS predictor
W = zeros(Nout, Nout);
Wo = [];
% compute first iteration with Xi(0) = Xi0
n = 1;
Wo = [Wo W];
xp(n, :) = Xi0 * W';
e(n, :) = Xi(n, :) - xp(n, :);
ne(n) = norm(e(n, :));
W = W + rp.mu * e(n, : )' * Xi0;
for n = 2:N,
% save W matrix
Wo = [Wo W];
% predict next sample and error
xp(n, :) = Xi(n-1, :) * W';
e(n, :) = Xi(n, :) - xp(n, :);
ne(n) = norm(e(n, :));
if (rp.verbose ~= 0)
disp(['time step ', int2str(n), ': mag. pred. err. = ', num2str(ne(n))]);
end;
% adapt weight matrix and step size
W = W + rp.mu * e(n, :)' * Xi(n-1, :);
if (rp.decay == 1)
rp.mu = rp.mu * n/(n+1); % use O(1/n) decay rate
end;
end