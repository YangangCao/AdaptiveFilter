clear all; clc; close all;

load 03_bbc_history_org.mat
clear sizeinfo
y=y(1:30*Fs);
h=firpm(160,[0 3.5 4.4 22.05]/22.05,[1 1 0 0]);
y=filter(h,1,y);
y=y(1:5:end);
Fs=Fs/5;
%%
Lh=80;
h=firpm(Lh,[0 0.47 0.53 1],[1 1 0 0]);
h=h.*(1j.^[0:Lh]);
y=filter(h,1,y);
y=2*real(y.*exp(1j*(0/Fs)*2*pi*[0:length(y)-1]'));
x=audioplayer(y,Fs);
play(x)
