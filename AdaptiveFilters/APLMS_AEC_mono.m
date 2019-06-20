%clear all; clc; close all;

load 01_memory
load wo
y=filter(wo,1,x);
eo=0.001*randn(size(y));
yeo=y+eo;
r=y;
%sound(y,Fs)
w=zeros(Lw,1);
e=y;
zeta=zeros(size(y));
mu=0.5;
M=5;
X=zeros(Lw,M);
for n=Lw:Lw+M-2
    xtdl=x(n:-1:n-Lw+1);
    X=[xtdl X(:,1:M-1)];
end
for n=Lw+M-1:length(x)-M
    xtdl=x(n:-1:n-Lw+1);
    X=[xtdl X(:,1:M-1)];
    Yhat=X'*w;
    yhat=Yhat(1);
    E=yeo(n:-1:n-M+1)-Yhat;
    e(n)=E(1);
    r(n)=y(n)-yhat;
    w=w+mu*X*inv(X'*X+0.1*eye(M))*E;
    zeta(n)=(w-wo)'*(w-wo);
end
E=audioplayer(e,Fs);
play(E)
ERLE=10*log10((filter(.001,[1 -.999],y.^2)+1e-12)./(filter(.001,[1 -.999],r.^2)+1e-12));
figure(2),plot(ERLE)
ERLEAPLMS=ERLE;
eAPLMS=e;
save ERLEAPLMS ERLEAPLMS eAPLMS
