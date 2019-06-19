function [Wo, xp, gamma, e] = qrd_rls_AR_pred(Xi, Y, verbose, lambda, init, init_val)
% qrd_rls_AR_pred.m
% function [Wo, xp, gamma, e] = qrd_rls_AR_pred(Xi, Y, verbose, lambda, init, init_val)
%
% qrd_rls_AR_pred.m - use the QR decomposition-based RLS algorithm
% to predict complex-valued AR process
% written for MATLAB 4.0
%
% Reference: Haykin, _Adaptive Filter Theory_, 2nd (corr.) ed., 1991
%
%
% Input Parameters
% Xi : matrix of training/test points - each row is
% considered a datum
% Y : vector of corresponding desired outputs for
% predictor
% verbose : set to 1 for interactive processing
% lambda : the initial value of the forgetting factor
% init,
% init_val: initialization method and value
% mvdr' - init_val is steering vector
%
%
% Output Parameters:
% Wo : row-wise matrix of Hermitian transposed weights
% at each iteration
% xp : row vector of predicted outputs
% gamma : row vector of a priori to a posteriori conversion factors
% e : row vector of a posteriori prediction errors Y - xp
Nout = size(Xi, 2);
% length of maximum number of timesteps that can be predicted
N = size(Xi, 1);
% order of predictor
M = size(Xi, 2);
% initialize weight matrix and associated parameters
% for QRD RLS predictor
W = zeros(Nout, 1);
Wo = [];
gamma = [];
lambda_root = sqrt(lambda);
Phi_root = zeros(M, M);
p_H = zeros(1, M);
W_H = zeros(1, M);
z = zeros(1, M);
% initialize predictor until Phi_root is full-rank '
for n = 1:M,
    % compute post-array from pre-array via QRD. Note that computation
    % of p_H, eta, and W_H are omitted in the initialization period.
    u = Xi(n, :).';
    pre_array = [ [lambda_root*Phi_root u] ; [lambda_root*p_H Y(n)] ; [z 1] ];
    post_array = triu(qr(pre_array'))';
    Phi_root = post_array(1:M, 1:M);
    gamma_root = post_array(M+2, M+1);
    % save results
    Wo = [ Wo; W_H ];
    gamma = [gamma gamma_root^2 ];
    % predict next sample and compute error
    xp(n) = W_H * u;
    e(n) = Y(n) - xp(n);
    if (verbose ~= 0)
        disp(['time step ', int2str(n), ': mag. pred. err. = ' , num2str(abs(e(n)))]);
    end;
end % for n
% if initializing for MVDR adaptive beamforming, set the
% auxiliary vector p after Phi_root is full-rank
if (strcmp(init, 'mvdr')),
    p_H = (Phi_root \ init_val)';
end;
% run normally after initialization is complete
for n = (M+1):N,
    % compute post-array from pre-array via QRD
    u = Xi(n, :).';
    pre_array = [ [lambda_root*Phi_root u] ; [lambda_root*p_H Y(n)] ; [z 1] ];
    post_array = triu(qr(pre_array'))';
    Phi_root = post_array(1:M, 1:M);
    p_H = post_array(M+1, 1:M);
    ganuna_root = post_array(M+2, M+1);
    eta = post_array(M+1, M+1) / gamma_root;
    W_H = (Phi_root' \ p_H')';
    if (strcmp(init, 'mvdr')), W_H = W_H / (p_H * p_H'); end;
    % save results
    Wo = [Wo; W_H];
    gamma = [gamma gamma_root^2 ];
    % predict next sample and compute error
    xp(n) = W_H * u;
    e(n) = Y(n) - xp(n);
    if (verbose ~= 0)
        disp(['time step ', int2str(n), ': mag. pred. err. = ', num2str(abs(e(n)))]);
    end;
end % for n