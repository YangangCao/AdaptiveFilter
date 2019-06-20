%
% Direction of arrival based on the MUSIC method
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
sigma_nu=sqrt(0.01);
X=[];
for n=1:100
    X=[X sqrt(P0)*randn*s0+sqrt(P1)*randn*s1+sqrt(P2)*randn*s2...
        +sigma_nu*randn(M,1)];
end
R=X*X'/100;
[V,D]=eig(R);
Qnu=V(:,1:M-3);
Qs=V(:,M-2:M);
for theta=-50:50
    s=exp(1j*pi*[0:M-1]'*sin(theta*pi/180));
    vnu=Qnu'*s;
    vs=Qs'*s;
    S(theta+51)=1/(vnu'*vnu);
    %S(theta+51)=vs'*vs;
end
S=S/max(S);
figure,axes('position',[0.25 0.25 0.5 0.5])
plot([-50:50],S,'k'),hold on
plot(20*[1 1],[0 1],'--')
text(20.5,0.9,'\theta_0')
plot(25*[1 1],[0 1],'--')
text(25.5,0.9,'\theta_1')
plot(-30*[1 1],[0 1],'--')
text(-29.5,0.9,'\theta_2')
xlabel('\theta')
hold off
