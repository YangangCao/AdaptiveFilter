%
% This MATLAB script provides an skeleton for simulation of a quadratue
% amplitude modulated (QAM) transmission system.
%
%%% PARAMETERS %%%
T=0.0001;       % Symbol/baud period 
L=100;          % Number of samples per symbol period
Ts=T/L;         % A dense sampling interval (approx to continuous time)
fc=100000;      % Carrier frequency at the transmitter
Dfc=0;          % Carrier frequency offset 
                % (carrier at transmitter minus carrier at receiver)
phic=0;         % Carrier phase offset
alpha=0.5;      % Roll-off factor for the square-root raised cosine filter 
sigma_v=0;      % Standard deviation of channel noise
c=1;            % Channel impulse response
%%% 4QAM symbols %%%
N=1000; s=sign(randn(N,1))+1i*sign(randn(N,1));
%%% Transmitter filterring %%%
pT=sr_cos_p(6*L,L,alpha);       % Transmit filter
xbbT=conv(expander(s,L),pT);    % Band-limited baseband transmit signal
%%% MODULATION %%%
t=[0:length(xbbT)-1]'*Ts;       % dense set of time points
xT=2*real(exp(i*2*pi*fc*t).*xbbT);
%%% CHANNEL (including aditive noise) %%%
xR=conv(c,xT); xR=xR+sigma_v*randn(size(xR)); % Received signal
%%% DEMODULATION %%%
t=[0:length(xR)-1]'*Ts;         % dense set of time points
xbbR=exp(-i*(2*pi*(fc+Dfc)*t-phic)).*xR;
%%% Receiver filtering  %%%
pR=pT; x=conv(xbbR,pR);



