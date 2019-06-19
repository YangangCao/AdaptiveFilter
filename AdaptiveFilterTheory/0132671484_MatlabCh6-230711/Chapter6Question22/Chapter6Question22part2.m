clear all;
close all;

mu = 1e-10; % step-size parameter 
g = 1; % unity gain

p = 5; % number of sensors
Ninit = p; % nuber of smaples needed for initialization
Nsnaps = [20,50,100]; % number of snapshots or iteration

mean_v = 0; % white noise mean
var_v = 1; % white noise variance

SNRdB = 10; % target signal to noise ratio
INRdB = 40; % interference signal to noise ratio

sin_theta = [-0.05 0]; % location of signal and interference

numst = 1000; % resolution in spatial response
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
    v = v_tmp(1, :)+ 1j*v_tmp(2, :); % additive white noise of desired mean/variance
    Psi = 2*pi*rand; % uniform random phase on interference
    Xi(i, :) = sig_x + A(2)*exp(1j*[1:p]*phi(2) + Psi) + v;
end; 
d = g*Xi(:,1);
u = diag(Xi(:,1))*(ones(Ndata,1)* e.')-Xi(:,2:p);

[W,xp] = lms(u,d,mu);
Wo = g- W * conj(e);
W = [Wo W]; 

% now, generate test vectors to compute spatially sampled response
W_H = conj(W(Ndata, :)); % Hemitian transpose of last one
st = linspace(-1,1,numst); % sine(theta) space
est = exp(-1j*pi*[0:(p-1)]'*st); % steering matrix
% amplitude response
P(:,n) = 20*log10(abs(W_H*est).^2); 
end

p1 = plot(st,P(:,1),'k',st,P(:,2),'r',st,P(:,3),'b');
%axis([-1 1 min(P) max(P)]);
axis([-1 1 -100 20]);
xlabel('sin \theta');
ylabel('Amplitude response (dB)');
legend('Amplitude response 5 adaption','Amplitude response 10 adaption'...
       ,'Amplitude response 20 adaption');
title(['MVDR beamfomer, numits = ',num2str(Nsnaps),',\mu = ',num2str(mu),...
      ', SNR =',num2str(SNRdB),', INR =',num2str(INRdB),', \phi target=',...
      num2str(sin_theta(1))]);
