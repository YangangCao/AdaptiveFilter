%MCMC equalizer
%for BPSK modulation only

clear
clc

% define channel 
%Prokais A channel 
% h=[0.04, -0.05,0.07,-0.21,-0.5,0.72,0.36,0,0.21,0.03,0.07 ];
h=[0.5 1 0.2];

L=length(h);            % channel length
N=2048;                % data bits
EbN0=-4:4:12;              % Eb/N0 points
EP=length(EbN0);        % number of Eb/N0 points
mc_iter = 10;           % Number of MCMC iterations

BER=zeros(EP);     % data bit error rate (filled during simulation)

% loop over all Eb/N0 values
for ep=1:EP
    
    % initialize random number generator
    rand('state',1063712);randn('state',5625914);
    
    % set noise variance
    nvar=h*h'/2/10^(EbN0(ep)/10);
    
    % data bit generation
    b=(rand(1,N)>0.5);
%      b(2:2:end) = 1;
    
    %BPSK modulation
    s = 1-2*b;
%     s= [s, zeros(1, L-1)];
    N1 = length(s)+L-1;
    x=conv(s,h)+randn(1, N1)*sqrt(nvar/2);
    
    %prior LLR for bits
    lambda2e= zeros(1, N);
    
    llr = mcmc_eq2(x, h, lambda2e,mc_iter, nvar);
    
    b2 = (llr<0);
    
    %Calculate BER
    ber(ep) = sum(abs(b2-b))/length(b);
end

semilogy(EbN0, ber);





    
    
    