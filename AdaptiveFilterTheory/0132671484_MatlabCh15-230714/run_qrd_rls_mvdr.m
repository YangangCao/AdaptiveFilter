function run_qrd_rls_mvdr(rp)
% run_qrd_rls_mvdr.m
% Computer Experiment
% Section 14.5, Adaptive Filter Theory, 3rd edition
% MVDR adaptive beamforming
% modify path to suit
% path(path, '/home/yee/aft/qrd_rls');
Ninit = rp.p;
Ndata = Ninit + rp.Nsnaps;
seed = 1;
lambda = 1;
% enter mean and variance of complex-valued AWGN
rp.mean_v = 0;
rp.var_v = 1;
% A_1, phi_1 are target signal amplitude/elec. angle
% A_2, phi_2 are interference signal amplitude/elec. angle
% s is steering vector along elec. angle of look direction of interest
A_1 = sqrt (rp.var_v) * 10^ (rp.TNRdB/20);
phi_1 = pi*rp.sin_theta_1;
A_2 = sqrt(rp.var_v) * 10^ (rp.INRdB/20);
phi_2 = pi*rp.sin_theta_2;
s = exp(-j*[0:(rp.p-1)]'*phi_1);
% setup input/output sequences
for i = 1:Ndata,
    % setup random disturbances
    randn('seed', i);
    vr = sqrt(rp.var_v/2) * randn(1, rp.p) + rp.mean_v;
    vi = sqrt(rp.var_v/2) * randn(1, rp.p) + rp.mean_v;
    v = vr + j*vi;
    rand('seed', i);
    %Psi = 2*pi*rand(1); % For Q 11.11
    Psi=0.2; %For Q 11.12
    Xi(i, :) = A_1*exp(j*[1:rp.p]*phi_1) + A_2*exp(j*[1:rp.p]*phi_2 + Psi) + v;
end;
Y = zeros(1, Ndata);
% run bea/nformer for indicated number of snapshots
[Wo, xp, gamma, e] = qrd_rls_AR_pred(Xi, Y, rp.verbose, lambda, 'mvdr', s);
% test vectors for spatially sampled response
W_H = Wo(Ndata, :);
st = -1:0.025:1;
est = exp(j*pi*[0:(rp.p-1)]'*st);
W=Wo; % simple renaming so that the format jives with thtat expected by the
% plot_mvdr.m routine
eval(['save ' rp.name])