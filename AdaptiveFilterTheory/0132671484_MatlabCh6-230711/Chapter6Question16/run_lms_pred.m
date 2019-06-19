function run_lms_pred(rp)
% rp is a structure of run parameters with elements
% Niter, Ndata, mult, verbose, alpha, a, var_v
% rp is created with the program makerp.m
% Computer Experiment
% Section 9.6, Adaptive Filter Theory, 3rd edition
% First-order prediction
seed = 0:(rp.Nruns-1);
rp.decay= 0;
Npred = rp.Ndata;
E = zeros(Npred*10, rp.Nruns);
WX = zeros(Npred, rp.Nruns);
Xi0 = 0;
for iter = 1:rp.Nruns,
    randn('seed', seed(iter));
    Xi = filter(1, [1 rp.a], [Xi0 ; sqrt(rp.var_v)*randn(rp.mult*rp.Ndata, 1)]);
    disp(['run # ' num2str(iter)]);
    disp([' covariance of AR process = ' num2str(cov(Xi))]);%
    Xi = Xi(1:(1000));
    lms_AR_pred;
    E(:, iter) = e(1:1000,:);
    Wx(:, iter) = Wo';
end;