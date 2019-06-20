%
% Direction of arrival based on the Bartlett method
% 
M=10;
theta0=20*pi/180;
theta1=25*pi/180;
theta2=-30*pi/180;
s0=exp(1j*pi*[0:M-1]'*sin(theta0));
s1=exp(1j*pi*[0:M-1]'*sin(theta1));
s2=exp(1j*pi*[0:M-1]'*sin(theta2));
P0=1;
P1=1;
P2=1;
sigma_nu=0.01;
X=[];
for n=1:100
    X=[X sqrt(P0)*randn*s0+sqrt(P1)*randn*s1+sqrt(P2)*randn*s2...
        +sigma_nu*randn(M,1)];
end
R=X*X'/100;
for theta=-50:50
    s=exp(1j*pi*[0:M-1]'*sin(theta*pi/180));
    S(theta+51)=s'*R*s;
end
figure,axes('position',[0.25 0.25 0.5 0.5])
plot([-50:50],S,'k'),hold on
plot(20*[1 1],[0 200],'--')
text(20.5,180,'\theta_0')
plot(25*[1 1],[0 200],'--')
text(25.5,180,'\theta_1')
plot(-30*[1 1],[0 200],'--')
text(-29.5,180,'\theta_2')
xlabel('\theta')
hold off
