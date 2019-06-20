%
% This script evaluates the minimum mean-squared error (MMSE) of the 
% symbol-spaced and half symbol-spaced equalizers. The equivalent baseband 
% channel h(t) is obtained according to (17.13). The MMSE values are then 
% calculated, following the formulations in Section 17.4. 
%
clear all,close all
%% parameters %%
T=0.0001; L=100; Ts=T/L; fs=1/Ts; fc=100000; 
alpha=0.5; sigmas=1; sigmanuc=0.01; 
c=[0.5 zeros(1,60) 1 zeros(1,137) 0.3];     % More common channel
c=[1 zeros(1,67) 0.75 zeros(1,145) 0.4];  % Less common channel 
c=c/sqrt(c*c');
pT=sr_cos_p(16*L,L,alpha)'; 
pR=pT; 
%% Construction of the equivalent baseband channel %%
p=conv(pT,pR);
c=c.*exp(-j*2*pi*[0:length(c)-1]*Ts*fc);
h=conv(c,p);
pR=sqrt(L/2)*pR(1:L/2:end);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%  T spaced equalizer   %
%%%%%%%%%%%%%%%%%%%%%%%%%%
N=21;       % equalizer length = N = 21
for k=1:L
    h0=h(k:L:end);
    Delta=round((length(h(1:L:end))+N)/2);
    C=toeplitz([h0 zeros(1,N-1)],[h0(1) zeros(1,N-1)]);
    P0=toeplitz([pR(1:2:end) zeros(1,N-1)],[pR(1) zeros(1,N-1)]);
    P1=toeplitz([pR(2:2:end) zeros(1,N-1)],[pR(2) zeros(1,N-1)]);
    Q=[C; (sigmanuc/sigmas)*P0; (sigmanuc/sigmas)*P1];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-14*eye(N));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmse0(k)=sigmas^2*real(1-w'*pyd);
    spower(k)=sum(abs(h0).^2);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%  T/2 spaced equalizer  %
%%%%%%%%%%%%%%%%%%%%%%%%%%
N=21;
hf=h(1:L/2:end);Delta=round((length(hf)+N)/4);
for k=1:L
    hf=h(k:L/2:end);
    C=toeplitz([hf zeros(1,N-1)],[hf(1) zeros(1,N-1)]); C=C(1:2:end,:);
    P=toeplitz([pR zeros(1,N-1)],[pR(1) zeros(1,N-1)]);
    Q=[C; (sigmanuc/sigmas)*P];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-6*eye(N));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmsef1(k)=sigmas^2*real(1-w'*pyd);
end

N=31;
hf=h(1:L/2:end);Delta=round((length(hf)+N)/4);
for k=1:L
    hf=h(k:L/2:end);
    C=toeplitz([hf zeros(1,N-1)],[hf(1) zeros(1,N-1)]); C=C(1:2:end,:);
    P=toeplitz([pR zeros(1,N-1)],[pR(1) zeros(1,N-1)]);
    Q=[C; (sigmanuc/sigmas)*P];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-6*eye(N));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmsef2(k)=sigmas^2*real(1-w'*pyd);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  T spaced decision fedback equalizer   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=21;       % FF filter length  
for k=1:L
    h0=h(k:L:end);
    Delta=round((length(h(1:L:end))+N)/2);
    C=toeplitz([h0 zeros(1,N-1)],[h0(1) zeros(1,N-1)]);
    M=length(h0)+N-1-Delta-1;  % FB filter length
    P0=toeplitz([pR(1:2:end) zeros(1,N-1)],[pR(1) zeros(1,N-1)]);
    P1=toeplitz([pR(2:2:end) zeros(1,N-1)],[pR(2) zeros(1,N-1)]);
    Q=[[C [zeros(Delta+1,M); eye(M)]]; ...
        [(sigmanuc/sigmas)*P0 zeros(length(P0(:,1)),M)]; ... 
        [(sigmanuc/sigmas)*P1 zeros(length(P1(:,1)),M)]];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-14*eye(N+M));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmsedf(k)=sigmas^2*real(1-w'*pyd);
end
%% Plots %%
figure(1)
subplot(3,2,1),plot([0:L-1]/L,spower); 
xlabel('Timing Phase'),ylabel('Signal Power')
subplot(3,2,3),semilogy([0:L-1]/L,mmse0,'k-',[0:L-1]/L,mmsef1,'k--',...
    [0:L-1]/L,mmsef2,'k-.',[0:L-1]/L,mmsedf,'k-')
xlabel('Timing Phase'),ylabel('MMSE'),
legend('T spaced equalizer (N=21)','T/2 spaced equalizer (N=21)',...
    'T/2 spaced equalizer (N=31)','DF equalizer (N=21)','Location','BestOutside')
[mx,k]=max(mmse0);
h0max=h(k:L:end);H0max=abs(fft(h0max,1000)); F=[0:999]/1000;
[mx,k]=min(mmse0);
h0min=h(k:L:end);H0min=abs(fft(h0min,1000));
subplot(3,2,5),plot(F,H0min,'k-',F,H0max,'k--')
legend('optimum \tau','worst \tau','Location','BestOutside')
xlabel('Frequency (\omega/2\pi)')
ylabel('Magnitude')
