close all
N=32;

Df=4/N;
W=Df/2;
n=[1:N-1];
r0=[2*W sin(2*pi*W*n)./(pi*n)];
R=toeplitz(r0);
[v,lambda]=eig(R);
v=v(:,end);v=v/sum(v);
V=20*log10(abs(fft(v,501)));
V4=fftshift(V);
n=[0:N-1];
wHanning=0.54-0.46*cos(2*pi*n/(N-1));wHanning=wHanning/sum(wHanning);
WHanning=fftshift(20*log10(abs(fft(wHanning,501))));
f=[-250:250]/500;
axes('position',[0.1 0.3 0.35 0.35])
plot(f,WHanning,'--k',f,V4,'-k')
legend('Hamming','prolate')%,'Location','SouthEast')
axis([-0.5 0.5 -120 20])
xlabel('NORMALIZED FREQUENCY, f')
ylabel('MAGNITUDE, dB')
hold on

Df=6/N;
W=Df/2;
n=[1:N-1];
r0=[2*W sin(2*pi*W*n)./(pi*n)];
R=toeplitz(r0);
[v,lambda]=eig(R);
v=v(:,end);v=v/sum(v);
V=20*log10(abs(fft(v,501)));
V6=fftshift(V);
n=[0:N-1];
wBlackman=0.42-0.5*cos(2*pi*n/(N-1))+0.08*cos(4*pi*n/(N-1));
wBlackman=wBlackman/sum(wBlackman);
WBlackman=fftshift(20*log10(abs(fft(wBlackman,501))));
axes('position',[0.6 0.3 0.35 0.35])
plot(f,WBlackman,'--k',f,V6,'-k')
legend('Blackman','prolate')%,'Location','SouthEast')
axis([-0.5 0.5 -120 20])
xlabel('NORMALIZED FREQUENCY, f')
ylabel('MAGNITUDE, dB')

% axes('position',[0.25 0.25 0.5 0.5])
% plot(f,Wrect,'--k',f,V2,'-k',f,WHanning,'--k',f,V4,'-k',f,WBlackman,'--k',f,V6,'-k')
% legend('rectangular','prolate (\Deltaf=2/N)','Hanning','prolate (\Deltaf=4/N)','Blackman','prolate (\Deltaf=6/N')
% axis([-0.5 0.5 -120 20])
% xlabel('NORMALIZED FREQUENCY')
% ylabel('MAGNITUDE, dB')
% hold
% axes('position',[0.5 0.5 0.3 0.3])
% plot(f,V)
% axis([0 4*W/16 -80 20])

