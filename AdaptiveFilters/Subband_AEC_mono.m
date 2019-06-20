%
% Subband AEC
%
%%
% parameters
K_a=5;
K_s=3;
M_a=19;
M_s=32;M=M_s;
Delta_a=K_a*M_a;
Delta_s=K_s*M_s;
N_a=2*K_a*M_a+1;
N_s=2*K_s*M_s+1;
Delta=ceil((Delta_a+Delta_s)/M)*M;
alpha_a=3/16;
alpha_s=7/19;
L=16;
%%
% filter design
h_a=eigenfir(M_a,N_a,alpha_a,Delta_a)';
h_s=eigenfir(M_s,N_s,alpha_s,Delta_s)';
%%
load 01_memory
load wo
y=filter(wo,1,x);
eo=0.001*randn(size(y));
yeo=y+eo;
%sound(y,Fs)
Lw_sub=ceil(Lw/L);
w=zeros(Lw_sub,M_s);
e=0;
r=0;
yhat=0;
zeta=zeros(size(y));
mu=0.5;
%%
for i=0:M/2
    h_a_i=h_a.*exp(1j*(2*pi*i/M)*[0:length(h_a)-1]');
    h_s_i=h_s.*exp(1j*(2*pi*i/M)*[0:length(h_s)-1]');
    x_sub=decimator(filter(h_a_i,1,x),L);
    yeo_sub=decimator(filter(h_a_i,1,yeo),L);
    y_sub=decimator(filter(h_a_i,1,y),L);
    w_sub=zeros(Lw_sub,1);
    e_sub=y_sub;
    r_sub=y_sub;
    y_sub_hat=y_sub;
    for n=Lw_sub:length(x_sub)
        xtdl=x_sub(n:-1:n-Lw_sub+1);
        y_sub_hat(n)=w_sub'*xtdl;
        e_sub(n)=yeo_sub(n)-y_sub_hat(n);
        r_sub(n)=y_sub(n)-y_sub_hat(n);
        w_sub=w_sub+mu/(xtdl'*xtdl+0.00001)*xtdl*e_sub(n)';
    end
    e=e+real(filter(h_s_i,1,expander(e_sub,L)));
    r=r+real(filter(h_s_i,1,expander(r_sub,L)));
    yhat=yhat+real(filter(h_s_i,1,expander(y_sub,L)));
end
%%
e=2*L*e;
yhat=2*L*real(yhat);
r=2*L*r;
% e=[zeros(Delta,1); yeo(1:end-Delta)]-yhat(1:length(yeo));
% r=[zeros(Delta,1); y(1:end-Delta)]-yhat(1:length(y));
E=audioplayer(e,Fs);
play(E)
figure(1),plot(e)
ERLE=10*log10((filter(.001,[1 -.999],yhat.^2)+1e-12)./(filter(.001,[1 -.999],r.^2)+1e-12));
figure(3),plot(ERLE)
ERLESB=ERLE;
eSB=e;
save ERLESB ERLESB eSB
