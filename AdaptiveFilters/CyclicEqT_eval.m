%
% This program evaluates the minimum mean-squared error (MMSE) of the symbol-spaced 
% cyclic equalizer, according to the derivations in Section 11.5.1. 
% The equivalent baseband channel cBB(t) is obtained according to (3.54). 
% The MMSE values are then calculated, by convolving cBB and w.
%
% Both direct and indirect method of calculating equalizer tap weights are
% presented.
%
% clear all
Tb=0.0001; L=100; Ts=Tb/L; fs=1/Ts; fc=100000; 
N=16*L; alpha=0.5; sigmas=1; sigmanuc=0.01; 
TmgPhase=0.9; epsilon=1e-6;
c=[0.5 zeros(1,60) 1 zeros(1,137) 0.3];     % More common channel
c=[1 zeros(1,67) 0.75 zeros(1,145) 0.4];  % Less common channel 
c=c/sqrt(c*c');
pT=sr_cos_p(N,L,alpha)'; 
pR=pT; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Construction of the equivalent baseband channel using (3.54)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=conv(pT,pR);
c=c.*exp(-j*2*pi*[0:length(c)-1]*Ts*fc);
cBB=conv(c,p); cBB=cBB(round(TmgPhase*L+1):L:end);
pR=sqrt(L/2)*pR(1:L/2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Construction of s[n]       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=31;           % Equalizer order
N=N+1;          % Equalizer length
pilot=CycPilot(N-1);
s=pilot;
for k=1:3
    s=[s; s];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Construction of the matrix Y and direct evaluation of w    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=conv(s,cBB);
v=sigmanuc*(randn(2*length(y),1)+i*randn(2*length(y),1))/sqrt(2);
v=conv(pR,v);
y=y+v(1:2:2*length(y));
y0=y(160:-1:160-N+1);
Y=y0;
for k=1:N-1
    y0=[y0(end); y0(1:end-1)];
    Y=[Y y0];
end
w_direct=(Y*Y'+epsilon*eye(N))\(Y*conj(pilot));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find an estimate of channel response and sigmanuc  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=pilot;
s0=flipud(b);
S=s0;
for k=1:N-1
    s0=[s0(end); s0(1:end-1)];
    S=[S s0];
end
y0=(y(130:130+N-1)+y(130+N:130+2*N-1))/2;
cBBest=((1/(b'*b))*S*conj(y0))';
[cmax,indx]=max(abs(cBBest));
while indx~=round(N/2)
    cBBest=[cBBest(end) cBBest(1:end-1)];
    [cmax,indx]=max(abs(cBBest));
end
v=y(130:130+N-1)-y(130+N:130+2*N-1);
sigmanucest=sqrt(sum(abs(v).^2)/(2*N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Optimum equalizer, MMSE, and the MSE arising  %
%  form the designed cyclic equalizers arising   %
%  from direct and inderect methods.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta=round((length(cBB)+N)/2);
C=toeplitz([cBB zeros(1,N-1)],[cBB(1) zeros(1,N-1)]);
P0=toeplitz([pR(1:2:end) zeros(1,N-1)],[pR(1) zeros(1,N-1)]);
P1=toeplitz([pR(2:2:end) zeros(1,N-1)],[pR(2) zeros(1,N-1)]);
Q=[C; (sigmanuc/sigmas)*P0; (sigmanuc/sigmas)*P1];
d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
Ryy=(Q.'*conj(Q)+1e-12*eye(N));
pyd=(Q.'*conj(d));
wo=Ryy\pyd;                     % Optimum equalizer
mmse=sigmas^2*real(1-wo'*pyd);  % MMSE

% Here, MSE arising from direct method is calculated.
% Since y[n] and pilot may not be aligned, the following for loop
% circularly rotates w find the best alignment; the one that results in
% smallest value of MSE.
w=w_direct;
for k=1:N
    mse(k)=sigmas^2*real(w'*Ryy*w-2*real(w'*pyd)+d'*d);
    w=[w(end); w(1:end-1)];
end
mse_direct=min(mse);

cBB=cBBest;
sigmanuc=sigmanucest;
Delta=round((length(cBB)+N)/2);
C=toeplitz([cBB zeros(1,N-1)],[cBB(1) zeros(1,N-1)]);
P0=toeplitz([pR(1:2:end) zeros(1,N-1)],[pR(1) zeros(1,N-1)]);
P1=toeplitz([pR(2:2:end) zeros(1,N-1)],[pR(2) zeros(1,N-1)]);
Q=[C; (sigmanuc/sigmas)*P0; (sigmanuc/sigmas)*P1];
d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
Ryy1=(Q.'*conj(Q)+1e-12*eye(N));
pyd1=(Q.'*conj(d));
w=Ryy1\pyd1;

for k=1:N
    mse(k)=sigmas^2*real(w'*Ryy*w-2*real(w'*pyd)+d'*d);
    w=[w(end); w(1:end-1)];
end
mse_indirect=min(mse);

% [mmse mse_direct mse_indirect]