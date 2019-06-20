%This program implement the LMS ALE, using Algorithm 2
%of Chapter 10 (Table 10.3).
%
%
%All parameters are specifically set to obtain Fig. 10.9b.
%The parameters may be changed to simulate other scenarios
%by editing the program.
%
%The variables names follow those in the text.
%The letter "v" has been added to indicate vectors.
%
%
% Last updated on April 28, 1998
%

itn=15000;
u=zeros(size(1:itn));
y=[0 0 0];
alpha=y;
beta=y;
e=y;
s=0.25;
theta=pi/2;
mus=0.001;
muthta=0.1;
sv2=zeros(size(1:itn));
wv2=sv2;
thetav2=sv2;
thetao=pi/3;
dthetao=0;
thetaok=0;
for k=2:itn
	if k==2000
		dthetao=pi/15000
	end
	if k==6000
		dthetao=0
		thetao=pi/10
	end
	if k==10000
		thetao=pi/20
	end
	thetao=thetao+dthetao;
	thetaok=thetaok+thetao;
	if thetaok>2*pi
		thetaok=thetaok-2*pi;
	end
	w=cos(theta);
	u(k+1)=sin(thetaok)+0.707*randn; %rndu(k);
	alpha(2:3)=alpha(1:2);
	beta(2:3)=beta(1:2);
	e(2:3)=e(1:2);
	y(2:3)=y(1:2);
	y(1)=(1+s)*w*y(2)-s*y(3)+(1-s)*(w*u(k)-u(k-1));
	e(1)=u(k+1)-y(1);
	alpha(1)=(1+s)*w*alpha(2)-s*alpha(3)-sin(theta)*((1+s)*y(2)+(1-s)*u(k));
	theta=theta+muthta*(1-s)^3*e(1)*alpha(1);
	if theta<0.05
		theta=0.05;
	end
	beta(1)=(1+s)*w*beta(2)-s*beta(3)-(w*e(2)-e(3));
	s=s+(mus/(1-s)^2)*(y(1)*y(1)+(1-s*s)*y(1)*beta(1));
	if s<0.25
		s=0.25;
	elseif s>0.9
		s=0.9;
	end
	sv2(k)=s;
	wv2(k)=w;
	thetav2(k)=thetao;
end
nn=[1:5:itn];
plot(nn,sv2(nn),'-.',nn,acos(wv2(nn)),'--',nn,thetav2(nn),'-')

