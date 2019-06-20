% %clear all; clc; close all;
% 
% [x, Fs] = wavread('02_man_npr.wav');
% h=firpm(160,[0 3.5 4.4 22.05]/22.05,[1 1 0 0]);
% x=filter(h,1,x);
% x=x(:,1);x=x(1:5:end);
% Fs=Fs/5; x=x(1:Fs*10);
% x=[x; x; x];
% Playing the wav file
%sound(x, Fs);
%%
load 01_memory
%%
N=1024; M=8;
% wo=randn(N,1)./(50+[1:N]');
% wo=wo/sqrt(wo'*wo);  
load wo
Lw=length(wo);
y=filter(wo,1,x);
eo=0.001*randn(size(y));
yeo=y+eo;
r=y;
e=y;
%%
P=ones(M+1,1);
w=zeros(N,1);
kappa=zeros(M,1);
b=x(N+M:-1:N);
bp=zeros(M,1);           %'bp' denotes the vector b' of Table 11.7
%Similarly 'fp' is used below for f'.
u_a=zeros(N,1);
beta=0.98;
alpha=0.9;
mu_po=0.001;
mu=0.5;
epsilon=0.05;
for n=N+M+1:length(x)
    %
    %	Lattice Predictor
    %
    b_old=b;
    f(1)=x(n);
    b(1)=f(1);
    P(1)=beta*P(1)+0.5*(1-beta)*(f(1)^2+b_old(1)^2);
    for m=1:M
        f(m+1)=f(m)-kappa(m)*b_old(m);
        b(m+1)=b_old(m)-kappa(m)*f(m);
        kappa(m)=kappa(m)+2*mu_po*(P(m)+epsilon)^(-1)*(f(m)*b(m+1)+b_old(m)*f(m+1));
        P(m+1)=beta*P(m+1)+0.5*(1-beta)*(f(m+1)^2+b_old(m+1)^2);
        if abs(kappa(m))>alpha
            kappa(m)=sign(kappa(m))*alpha;
        end
    end
    %
    %  u_a(n) update
    %
    u_a(2:N)=u_a(1:N-1);
    b_old=bp;
    fp(1)=b(M+1);
    bp(1)=fp(1);
    for m=1:M
        fp(m+1)=fp(m)-kappa(m)*b_old(m);
        bp(m+1)=b_old(m)-kappa(m)*fp(m);
    end
    u_a(1)=fp(M+1)/(P(M+1)+epsilon);
    %
    %   Filtering
    %
    xtdl=x(n-M:-1:n-N-M+1);
    yhat=w'*xtdl;
    e(n)=yeo(n-M)-yhat;
    w=w+(mu/(abs(u_a'*xtdl)+0.1))*e(n)*u_a;%
    r(n)=y(n-M)-yhat;
    %xi(n)=xi(n)+e*e;
end
%%
E=audioplayer(e,Fs);
play(E)
figure(1), plot(e)
%figure(2),plot(y)
ERLE=10*log10((filter(.001,[1 -.999],y.^2)+1e-12)./(filter(.001,[1 -.999],r.^2)+1e-12));
t=[0:length(ERLE)-1]/Fs;
figure(3),plot(t,ERLE)
ERLELMSNewton=ERLE;eLMSNewton=e;
save ERLELMSNewton ERLELMSNewton eLMSNewton

