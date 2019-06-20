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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Early-late gate timing recovery %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0;             % Should be changed for the modified algorithm             
mu0=0.005;          % Step-size parameter
dtau=12;            % (delta tau) times L
mu=mu0*(L/4)/dtau;  % Adjusted step-size
kk=1;xp=0;xm=0;
start=5*L+1;        % To drop the transient of x(t) at the beginning
tau=0.3*ones(1,floor((length(x)-start)/L)); % Initialize the timing offset
            % The timing offset tau is adjusted as the algorithm proceeds.
for k=start:L:length(tau)*L
    tauT=round(tau(kk)*L);
    xp=sqrt(1-beta^2)*x(k+tauT+dtau)-beta*xp;
    xm=sqrt(1-beta^2)*x(k+tauT-dtau)-beta*xm;
    tau(kk+1)=tau(kk)+mu*(abs(xp)^2-abs(xm)^2);
    kk=kk+1;
end
figure, axes('position',[0.1 0.25 0.7 0.5]), plot(tau(1:kk),'k')
xlabel('Iteration Number, n'), ylabel('\tau(n)')


