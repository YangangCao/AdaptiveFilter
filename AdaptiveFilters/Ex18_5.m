%
% Example 18.5
%
M=10;
P0=4;
P1=1;
P2=1;
a=0.95;
sigma_delta=2;
sigmanu2=0.01;
theta0=45*pi/180;
theta1=75*pi/180;
theta2=-30*pi/180;
so=exp(1i*pi*sin(theta0)*[0:M-1]');
s0=exp(1i*a*pi*sin(theta0)*[0:M-1]');
s1=exp(1i*a*pi*sin(theta1)*[0:M-1]');
s2=exp(1i*a*pi*sin(theta2)*[0:M-1]');
R=P0*s0*s0'+P1*s1*s1'+P2*s2*s2'+sigmanu2*eye(M);
p=0;
for k= 1:10000
    theta=(45+sigma_delta*randn)*pi/180;
    so=exp(1i*pi*sin(theta)*[0:M-1]');
    R=R+so*so';
    p=p+so;
end
woc=R\p;
for k=0:360
    theta=k*pi/180;
    G(k+1)=woc'*exp(1j*a*pi*[0:M-1]'*sin(theta));
end
theta=[0:360]*pi/180;
figure(1),polar(theta,abs(G))
figure(2),plot(theta*180/pi,abs(G))