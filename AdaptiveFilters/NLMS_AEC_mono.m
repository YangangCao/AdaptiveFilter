%clear all; clc; close all;

load 01_memory
load wo
Lw=length(wo);
figure(1),plot(wo,'r')
y=filter(wo,1,x);
eo=0.001*randn(size(y));
d=y+eo;
r=y;
%sound(y,Fs)
w=zeros(Lw,1);
e=y;
zeta=zeros(size(y));
mu=0.5;psi=0.1;alpha=0.995;
Pd=1; Pyhat=1; Pe=1;
flag=1; % flag=0 implies fixed step-size
        % flag=1 implies variable ste-size
for n=length(w):length(x)
    xtdl=x(n:-1:n-Lw+1);
    yhat=w'*xtdl;
    e(n)=d(n)-yhat;
    r(n)=y(n)-yhat;
    if flag==1
        Pd=alpha*Pd+(1-alpha)*d(n)*d(n);
        Pyhat=alpha*Pyhat+(1-alpha)*yhat*yhat;
        Pe=alpha*Pe+(1-alpha)*e(n)*e(n);
        mu=1-0.5*Pyhat/Pd;
        if mu>1
            m=1;
        elseif mu<0
            mu=0;
        end
    end
    w=w+mu/(xtdl'*xtdl+psi)*xtdl*e(n);
    zeta(n)=(w-wo)'*(w-wo);
end
E=audioplayer(e,Fs);
play(E)
figure(1),plot(e)
ERLE=10*log10((filter(.001,[1 -.999],y.^2)+1e-12)./(filter(.001,[1 -.999],r.^2)+1e-12));
t=[0:length(ERLE)-1]/Fs;
figure(4),plot(t,ERLE)
%figure(5),semilogy(zeta)
ERLENLMS=ERLE;
eNLMS=e;
save ERLENLMS ERLENLMS eNLMS
