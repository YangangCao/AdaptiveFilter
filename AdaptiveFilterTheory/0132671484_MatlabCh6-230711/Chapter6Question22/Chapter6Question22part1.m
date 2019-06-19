clc
clear

SignalToNoise=20;%dB
InterferenceToNoise=40;%dB
mu=1e-5;
var=1;
Ndata=100+5;
count=[1:5]-1;
target=asin(-0.2);
Interference=asin(0);
A_1mag=sqrt(var)*10^(SignalToNoise/20); %signal to noise as a ratio
phi_1=pi*asin(0);
A_2mag=sqrt(var)*10^(InterferenceToNoise/20);   %interferencenoise as ratio
phi_2=pi*asin(-0.2);
s=exp(1i*count*phi_2); %steering vector

for n=1:Ndata
    vr=sqrt(var/2)*randn(1,5);
    vi=sqrt(var/2)*randn(1,5);
    v=vr+1i*vi;
    Psi=2*pi*rand';
    A_2=A_2mag*randn(1);
    A_1(n)=A_1mag*randn(1);
    Xi(n,:)=A_1(n)*exp(1i*phi_1*count)+A_2*exp(1i*(count)*phi_2)+v;
end

d=A_1;
%u=diag(Xi(:,1))*(ones(Ndata,1)*e)-Xi(:,2:5);

Weights=zeros(5,1);
for ell=1:Ndata
    error(ell)=d(ell)-Xi(ell,:)*Weights;
    error(ell)=norm(error(ell))/d(ell);
    Weights=Weights+mu*error(ell)*Xi(ell,:)';
end
k=-pi:0.01:pi;
theta=pi*sin(k);
plot(k,10*log10(real([ones(size(k,2),1),exp(-1j*theta'),exp(-2j*theta'),exp(-3j*theta'),exp(-4j*theta')])*Weights))
%plot(error); %steering vector)
