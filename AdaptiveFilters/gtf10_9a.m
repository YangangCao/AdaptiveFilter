%This program implement the LMS ALE, using Algorithm 1
%of Chapter 10 (Table 10.2).
%
%All parameters are specifically set to obtain Fig. 10.9a.
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
w=0;
mus=0.001;
muw=0.05;
sv1=zeros(size(1:itn));
wv1=sv1;
thetav1=sv1;
theta=pi/3;
dtheta=0;
thetak=0;
for k=2:itn
	if k==2000
		dtheta=pi/15000
	end
	if k==6000
		dtheta=0
		theta=pi/10
	end
	if k==10000
		theta=pi/20
	end
	theta=theta+dtheta;
	thetak=thetak+theta;
	if thetak>2*pi
		thetak=thetak-2*pi;
	end
	u(k+1)=sin(thetak)+0.707*randn;	%rndu(k);
	alpha(2:3)=alpha(1:2);
	beta(2:3)=beta(1:2);
	e(2:3)=e(1:2);
	y(2:3)=y(1:2);
	y(1)=(1+s)*w*y(2)-s*y(3)+(1-s)*(w*u(k)-u(k-1));
	e(1)=u(k+1)-y(1);
	alpha(1)=(1+s)*w*alpha(2)-s*alpha(3)+(1+s)*y(2)+(1-s)*u(k);
	w=w+muw*(1-s)^3*e(1)*alpha(1);
	if abs(w)>0.999
		w=0.999*sign(w);
	end
	beta(1)=(1+s)*w*beta(2)-s*beta(3)-(w*e(2)-e(3));
	s=s+(mus/(1-s)^2)*(y(1)*y(1)+(1-s*s)*y(1)*beta(1));
	if s<0.25
		s=0.25;
	elseif s>0.9
		s=0.9;
	end
	sv1(k)=s;
	wv1(k)=w;
	thetav1(k)=theta;
end
nn=[1:5:itn];
plot(nn,sv1(nn),'-.',nn,acos(wv1(nn)),'--',nn,thetav1(nn),'-')


