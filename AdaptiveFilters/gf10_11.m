%This program implement the LMS ALE, using Algorithm 2
%of Chapter 10 (Table 10.3).
%The input is the summation of a number of sine-waves, 
%plus white noise.
%
%
%All parameters are specifically set to obtain Fig. 10.11.
%The parameters may be changed to simulate other scenarios
%by editing the program.
%
%To run this program on the student edition of MATLAB, reduce 
%the iteration number "itn" to lower value, such as itn=3200.
%
%
% Last updated on April 28, 1998
%

itn=10000;       %No. of iterations.
M=4;		        %No. of spectral lines.
u=zeros(M+1,3);
y=u;
alpha=u;
beta=u;
s=0.25*ones(1,M);
theta=(pi/2)*ones(1,M);
mus=0.001*ones(1,M);
muthta=0.02*ones(1,M);
Pu=ones(1,M);
sm=zeros(M,itn);
wm=sm;
thetao=[pi/1.8 pi/3.5 pi/6 pi/12];
phases=[0 pi/7 pi pi/3.76];
amps=[1 2 0.25 0.5];
um=zeros(M+1,itn);
for k=1:itn
   w=cos(theta);
   u(:,2:3)=u(:,1:2);
   y(:,2:3)=y(:,1:2);
   alpha(:,2:3)=alpha(:,1:2);
   beta(:,2:3)=beta(:,1:2);
   u(1,1)=sum(amps.*sin(thetao*k+phases))+0.5*randn;
   w=cos(theta);
   flag=ones(1,M);
   for l=1:M
	Pu(l)=0.02*u(l,1)*u(l,1)+0.98*Pu(l);
	y(l,1)=(1+s(l))*w(l)*y(l,2)-s(l)*y(l,3)+(1-s(l))*(w(l)*u(l,2)-u(l,3));
	u(l+1,1)=u(l,1)-y(l,1);
	alpha(l,1)=(1+s(l))*w(l)*alpha(l,2)-s(l)*alpha(l,3) ...
			-sin(theta(l))*((1+s(l))*y(l,2)+(1-s(l))*u(l,2));
	theta(l)=theta(l)+flag(l)*(muthta(l)/Pu(l))*((1-s(l))^3)*u(l+1,1)*alpha(l,1);
	if theta(l)<0
		theta(l)=0.1;
	end
	beta(l,1)=(1+s(l))*w(l)*beta(l,2)-s(l)*beta(l,3) ...
						-(w(l)*u(l+1,2)-u(l+1,3));
	s(l)=s(l)+flag(l)*(mus(l)/(Pu(l)*(1-s(l))^2))*(y(l,1)*y(l,1)+(1-s(l)*s(l))*y(l,1)*beta(l,1));
	if s(l)<0.25
		s(l)=0.25;
	elseif s(l)>0.9
		s(l)=0.9;
	end
	if s(l)<0.85
		flag(l+1:M)=zeros(size(l+1:M));
	end
   end
   sm(:,k)=s';
   wm(:,k)=w';
   um(:,k)=u(:,1);
end
figure(1),plot(wm')
figure(2),plot(sm')
