%
% OFDM transceiver
%
clear all, close all
%%%%%%%%%%%%%%%%
%  Parameters  %
%%%%%%%%%%%%%%%%
N=64;           % length of FFT/IFFT
Ncp=16;         % length of CP
Nactive=52;     % Number of active subcarriers
Nsymbols=50;    % Number of OFDM symbols, to transmit
Tb=0.0001;      % Symbol/baud period 
L=40;           % Number of samples per symbol period/oversampling factor
Ts=Tb/L;        % Sampling period
fc=50000;       % Carrier frequency at the transmitter
Dfc=0;          % Carrier frequency offset
phic=0;         % Carrier phase offset
K=8;            % Length of transmit and receive filters in Tb. 
sigmav=0.0001;       % Standard deviation of channel noise
c=1;   % Channel impulse response
%c=[1; zeros(47,1); 0.5; zeros(28,1); 0.8];
%%%%%%%%%%
% Source %
%%%%%%%%%%
b=sign(randn(2*Nactive*Nsymbols,1));
%%%%%%%%%
% coder %
%%%%%%%%%
sfreq=(b(1:2:end)+i*b(2:2:end))/sqrt(2);   % convert bits to QPSK symbols
                                 % These are in the frequency domain.
stime=OFDMmod(sfreq,N,Nactive,Ncp); % This functions duilds the OFDM symbols,
                                    % performs subcarriers modulation (IFFT),
                                    % and add cyclic prefix.

%%%%%%%%%%%%%%%%%%%
%   Add preamble  %
%%%%%%%%%%%%%%%%%%%
sp=0.25*sqrt(13/6)*(1+1i)*ifft([0 -1 -1 1 1 1 1 0 0 0 1 -1 1 -1 -1 1 ]');
shortpreamble=[sp; sp; sp; sp; sp; sp; sp; sp; sp; sp];
sl=ifft([ 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 ...
    0 0 0 0 0 0 0 0 0 0 0 ...
    1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1]');
longpreamble=[sl(end-31:end); sl; sl];
stime=[shortpreamble; longpreamble; stime];


%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
pT=firpm(K*L,[0 26 38 32*L]/(32*L),[1 1 0 0],[1 10])'; 
xbbT=conv(expander(stime,L),pT);

%%%%%%%%%%%%%%%%%%%%%%
%     MODULATION     %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xbbT)-1]'*Ts;      % Set the time indices
xT=real(exp(i*2*pi*fc*t).*xbbT);

%%%%%%%%%%%%%%%%%%%%%%
%      CHANNEL       %
%%%%%%%%%%%%%%%%%%%%%%
xR=conv(c,xT);
xR=xR+sigmav*randn(size(xR)); % Received signal

%%%%%%%%%%%%%%%%%%%%%%
%    DEMODULATION    %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xR)-1]'*Ts;     % Set the time indices
xbbR=2*exp(-i*(2*pi*(fc+Dfc)*t-phic)).*xR;

%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE FILTERING  %
%%%%%%%%%%%%%%%%%%%%%%
pR=pT;    
y=conv(xbbR,pR);
y=L*sqrt(N)*y(1:L:end);     % Decimation. 
        % The factor L*sqrt(N) is to bring up the signal
        % power to value around unity.

%%%%%%%%%%%%%%%%%%%%%%
% TIMING ACQUISITION %
%%%%%%%%%%%%%%%%%%%%%%
Ly=length(y);
for k=1:500
    ryy(k)=(y(k:k+N-1)'*y(k+N:k+2*N-1))/N;
end
figure,axes('position',[0.25 0.25 0.5 0.5])
plot(abs(ryy))
xlabel('n')
ylabel('r_{yy}(n)')
text(7,0.857,'n_1')
text(38,0.857,'n_2')
text(162,0.857,'n_3')
text(198,0.857,'n_4')

%%%%%%%%%%%%%%%%%%%%%%
% To be continued by the reader ...
%%%%%%%%%%%%%%%%%%%%%%