%
% This program generates a double talk situation and uses 
% power detection and coherence detection for double-talk 
% detection. The decision numbers D1 and D2, discussed in 
% the text are generated and presented in two separate
% figures.
%
%%
% Read a sample speech and reduce its sampling rate from 
% 44100 Hz to 44100/5 = 8820 Hz
%
load 01_memory_org
h=firpm(160,[0 3.5 4.4 22.05]/22.05,[1 1 0 0]);
y=filter(h,1,y);
x=y(1:5:end);
Fs=Fs/5; 
x=x(1:80*Fs);
%%
% Read a second speech signal, reduce its sampling rate,
% force to zero its sample values in the intervals (0 to 15)
% (30 to 45), (50 to 60) and (70 to 80) seconds. Then add this
% to the echo version of the first speech signal.
% This result in the double talk periods (15 to 30), (45 to 50),
% and (60 to 70) seconds.
% The duration of the speech signals is 80 seconds.
%
load 03_bbc_history_org
y=filter(h,1,y);
xd=y(1:5:end);
Fs=Fs/5;
xd(1:15*Fs)=zeros(size(xd(1:15*Fs)));
xd(30*Fs:45*Fs)=zeros(size(xd(30*Fs:45*Fs)));
xd(50*Fs:60*Fs)=zeros(size(xd(50*Fs:60*Fs)));
xd(70*Fs:80*Fs)=zeros(size(xd(70*Fs:80*Fs)));
xd=xd(1:80*Fs);
%
load wo
y=filter(wo,1,x);
d=y+xd;
%%
% Power based double-talk detector. 
% Produces the decision number D2.
% It is assumed that the echo canceler is perfectly adjusted,
% thus, y(n) is available.
% 
% Mean powers of y(n) and d(n) are canculated by average M=1000
% samples of the square of each.
%
M=1000;
sigmay2=filter(ones(M,1)/M,1,y.^2);
sigmad2=filter(ones(M,1)/M,1,d.^2);
D2=sigmay2./sigmad2;
%%
% Multitaper prototypes
%
L=256;K=8;
Df=1/L;
W=Df/2;
n=[1:L*K-1];
r0=[2*W sin(2*pi*W*n)./(pi*n)];
R=toeplitz(r0);
[v,lambda]=eig(R);
v=v(:,end-K+1:end); % the prototypes are in columns of v
%%
% Use Multitaper prototypes to extract narrow-band portions of
% x(n) and d(n). These will be then used to calculate the
% coherent function.
%
XK=zeros(L,floor(length(y)/L),K);
DK=XK;
for m=1:K
    D=d(1:L*floor(length(d)/L));
    D=reshape(D,L,length(D)/L); % arrange decimated samples in columns
    h=v(:,m);
    h=reshape(h,L,length(h)/L); % polyphase coefficients
    for k=1:L
        D(k,:)=filter(h(k,:),1,D(k,:));
    end
    DK(:,:,m)=fft(D);
    %
    X=x(1:L*floor(length(x)/L));
    X=reshape(X,L,length(X)/L); % arrange decimated samples in columns
    h=v(:,m);
    h=reshape(h,L,length(h)/L); % polyphase coefficients
    for k=1:L
        X(k,:)=filter(h(k,:),1,X(k,:));
    end
    XK(:,:,m)=fft(X);
end
%%
% Calculate the means across the data from different multitaper
% filter banks and across samples (here, 5 samples) in time.
%
DX=DK.*conj(XK);
DD=abs(DK).^2;
XX=abs(XK).^2;
DXmean=zeros(size(DX(:,:,1)));
DDmean=zeros(size(DX(:,:,1)));
XXmean=zeros(size(DX(:,:,1)));
for m=1:K
    DXmean=DXmean+DX(:,:,m);
    DDmean=DDmean+DD(:,:,m);
    XXmean=XXmean+XX(:,:,m);
end
for k=1:L
    DXmean(k,:)=filter(ones(1,5),1,DXmean(k,:));
    DDmean(k,:)=filter(ones(1,5),1,DDmean(k,:));
    XXmean(k,:)=filter(ones(1,5),1,XXmean(k,:));
end
%%
% The final step of the computation of the coherence function 
C=(abs(DXmean).^2)./(DDmean.*XXmean);
% Take the mean across frequency. This gives the decision number
% D1.
D1=mean(C(floor(L/50):floor(L/4),:));
%%
% Plot the results
%
t=[0:length(D1)-1]*L/Fs;
figure,axes('position',[0.25 0.25 0.5 0.5])
plot(t,D1)
xlabel('TIME (seconds)')
ylabel('D_1'), hold on
plot([15 15],[0 0.8],'--k')
plot([30 30],[0 0.8],'--k')
plot([45 45],[0 0.8],'--k')
plot([50 50],[0 0.8],'--k')
plot([60 60],[0 0.8],'--k')
plot([70 70],[0 0.8],'--k')
%
t=[0:length(D2)-1]/Fs;
figure,axes('position',[0.25 0.25 0.5 0.5])
plot(t,D2)
xlabel('TIME (seconds)')
ylabel('D_2'), hold on
plot([15 15],[0 1.4],'--k')
plot([30 30],[0 1.4],'--k')
plot([45 45],[0 1.4],'--k')
plot([50 50],[0 1.4],'--k')
plot([60 60],[0 1.4],'--k')
plot([70 70],[0 1.4],'--k')