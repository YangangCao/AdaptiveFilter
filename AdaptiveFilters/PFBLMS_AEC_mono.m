%clear all; clc; close all;

load 01_memory
% wo=randn(N,1)./(50+[1:N]');
% wo=wo/sqrt(wo'*wo);  
load wo
%%
N=1024; M=64; P=N/M;
Lw=length(wo);
y=filter(wo,1,x);
eo=0.001*randn(size(y));
yeo=y+eo;
r=y;
%%
wF=zeros(2*M,P);
xF=zeros(2*M,P);
mu=0.5;
e=y;
zeta=zeros(size(y));
for n=M+1:M:length(x)-M
    xF=[fft(x(n-M:n+M-1)) xF(:,[1:end-1])];
    yhat=ifft(sum((wF.*xF).').'); 
    yhat=real(yhat(M+1:end));
    E=yeo(n:n+M-1)-yhat;
    MU=mu*(sum((abs(xF).^2)')'+0.1).^(-1);
    EF=fft([zeros(M,1); E]);
    wF=wF+diag(MU.*EF)*conj(xF);
    waux=real(ifft(wF)); wF=fft([waux(1:M,:); zeros(M,P)]);
    e(n:n+M-1)=yeo(n:n+M-1)-yhat;
    r(n:n+M-1)=y(n:n+M-1)-yhat;
end
%%
E=audioplayer(e,Fs);
play(E)
ERLE=10*log10((filter(.001,[1 -.999],y.^2)+1e-12)./(filter(.001,[1 -.999],r.^2)+1e-12));
figure(3),plot(ERLE,'g')
ERLEFBLMS=ERLE;
eFBLMS=e;
save ERLEFBLMS ERLEFBLMS eFBLMS
