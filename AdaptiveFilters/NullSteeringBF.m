%
% Null steering beam former
%
clear all, close all
M=10;
theta0=20*pi/180;
theta1=25*pi/180;
theta2=-30*pi/180;
s0=exp(1j*pi*[0:M-1]'*sin(theta0));
s1=exp(1j*pi*[0:M-1]'*sin(theta1));
s2=exp(1j*pi*[0:M-1]'*sin(theta2));
S=[s0 s1 s2];
e0=[1 0 0]';
w=S*inv(S'*S)*e0;
for k=0:360
    theta=k*pi/180;
    G(k+1)=w'*exp(1j*pi*[0:M-1]'*sin(theta));
end
theta=[0:360]*pi/180;
polar(theta,abs(G))

