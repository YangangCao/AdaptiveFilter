clear all;
close all;
clc;
p = 5; % number of sensors
Ninit = p; % nuber of smaples needed for initialization
Nsnaps = [10, 10, 20]; % number of snapshots or iteration
mean_v = 0; % white noise mean
var_v = 1; % white noise variance
lambda = 0.999; % forgetting factor 
SNRdB = 10; % target signal to noise ratio
INRdB = 40; % interference signal to noise ratio
numst = 1000; % resolution in spatial response
sin_theta = [-0.05 0]; % location of signal and interference
phi = pi.*sin_theta; % equivalent electrical angle
A = sqrt(var_v)*10.^([SNRdB INRdB]./20); % parameter for target/interference signal amplitude, angle
% steering vector along electrical angle of look direction of interest
e = exp(-1j*[1:(p-1)]'*phi(1));
% setup input and output sequences

for n=1:3
    Ndata = Ninit + Nsnaps(n); % total number of data
    sig_x = A(1)*exp(1j*[1:p]*phi(1));
    for i = 1:Ndata,
    % random disturbances
        v_tmp = sqrt(var_v/2)*randn(2,p)+mean_v;
        v = v_tmp(1, :)+ 1j*v_tmp(2, :); % additive white noise 
        Psi = 2*pi*rand; % uniform random phase on interference
        Xi(i, :) = sig_x + A(2)*exp(1j*[1:p]*phi(2) + Psi) + v;
    end; 
    g = 1; % unity gain
    d = g*Xi(:,1);
    u = diag(Xi(:,1))*(ones(Ndata,1)* e.')-Xi(:,2:p);

    [W,xp] = rls(u,d,lambda);
    Wo = g- W * conj(e);
    W = [Wo W]; 

    % now, generate test vectors to compute spatially sampled response
    W_H = conj(W(Ndata, :)); % Hemitian transpose of last one
    st = linspace(-1,1,numst); % sine(theta) space
    est = exp(-1j*pi*[0:(p-1)]'*st); % steering matrix
    % amplitude response
    P(:,n) = 20*log10(abs(W_H*est).^2); 
end

plot(st,P(:,3),'k');
%axis([-1 1 min(P) max(P)]);
axis([-1 1 -100 40]);
xlabel('sin \theta');
ylabel('Amplitude response (dB)');
legend('Amplitude response 20 adaption');
title(['MVDR beamfomer,','\lambda = ',num2str(lambda),', SNR =',num2str(SNRdB),', INR =',num2str(INRdB),', \phi target=',num2str(sin_theta(1))]);
