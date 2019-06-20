%This is a VSLMS for estimation of a multipath channel.
%Here, we have considered a two path channel. The first
%path is fixed at \tau=2T. The second path starts at
%\tau=4T and slowly (over 100,000 symbols) moves (at constant rate) 
%to \tau=10T.
%
%Each path is assumed to be a slowly varying Gaussian process
%with an spectrum matching that of a single pole lowpass filter.
%The bandwidth of the latter filter is selected equal to the
%fade rate.
%
%Different fade rates may be set. For each particular
%run the two paths are assumed to fade with the same rate.
%
%
% Last updated on April 28, 1998
%

%
%       Set the channel specifications.
%
%npaths=2;
f_d=1;
f_s=2400;
deltao=0.02;
alpha=1-pi*f_d/f_s;
beta=sqrt(1-alpha^2);
itn=100000
N=16;
sigmaeo=sqrt(0.02);

muus=(0.1/N)*ones(N,1); 	%For 10 percent(approx.) misadjustment
gus2=ones(N,1);
wus=zeros(N,1);
gusold=wus;
xius=zeros(1,itn/100);
xiauxus=0;
fmuus=muus;
mu10us=zeros(1,itn/100);
	
musgn=(0.1/N)*ones(N,1); 	%For 10 percent(approx.) misadjustment
wsgn=zeros(N,1);
gsgnold=wsgn;
deltas=deltao;
xisgn=zeros(1,itn/100);
xiauxsgn=0;
%Musgn=zeros(N,itn/100);
fmusgn=musgn;
mu10sgn=zeros(1,itn/100);

muopt=(0.1/N)*ones(N,1); 	%For 10 percent(approx.) misadjustment
wopt=zeros(N,1);
gsgnold=wsgn;
xiopt=zeros(1,itn/100);
xiauxopt=0;
fmuopt=muopt;
mu10opt=zeros(1,itn/100);

a1=1;
a2=1;
tau1=2;
tau2=4;
dtau=10/itn;
t=[0:1:N-1]';
s=sign(randn(N,1));

for k=1:itn
%       
%   Update channel response.
%
	a1=beta*randn+alpha*a1;
	a2=beta*randn+alpha*a2;
	tau2=tau2+dtau;
	wo=a1*rcos50(t-tau1)+a2*rcos50(t-tau2);

%
%    Desired output
%
	d=wo'*s+sigmaeo*randn;

%
%    Un-sign (conventinal) VSLMS
%
	e=d-wus'*s;
	gus=e*s;
	wus=wus+2*muus.*gus;
	gus2=0.95*gus2+0.05*gus.*gus;
	z=deltao*gus.*gusold./gus2;
	for kk=1:N
		if abs(z(kk))>0.1
			z(kk)=0.1*sign(z(kk));
		end
	end
	muus=muus.*(1+z);
	if sum(muus)>1/3
		muus=muus/(3*sum(muus));
	end
	fmuus=0.01*muus+0.99*fmuus;
	gusold=gus;
	xiauxus=xiauxus+(wus-wo)'*(wus-wo);

%
%    Sign VSLMS
%
	e=d-wsgn'*s;
	gsgn=e*s;
	wsgn=wsgn+2*musgn.*gsgn;
	musgn=musgn.*(1+deltas*sign(gsgn).*gsgnold);
	if sum(musgn)>1/3
		musgn=musgn/(3*sum(musgn));
	end
	fmusgn=0.01*musgn+0.99*fmusgn;
	gsgnold=sign(gsgn);
	xiauxsgn=xiauxsgn+(wsgn-wo)'*(wsgn-wo);

%
%	LMS with optimum step-size parameters
%
	sigmaepo=beta*abs(rcos50(t-tau1)+rcos50(t-tau2));
	eta=(sum(sigmaepo)+sqrt(sum(sigmaepo)^2+4*sigmaeo^2))/2;
	muopt=sigmaepo/(2*eta);
	for kk=1:16
		if muopt(kk)<0.001
			musgn(kk)=0.001;
		end
	end
	if sum(muopt)>1/3
		muopt=muopt/(3*sum(muopt));
	end
	e=d-wopt'*s;
	wopt=wopt+2*e*muopt.*s;
	xiauxopt=xiauxopt+(wopt-wo)'*(wopt-wo);



	if k==1;
		xiopt(1)=wo'*wo;
		xius(1)=wo'*wo;
		xisgn(1)=wo'*wo;
	end
		
	if rem(k,100)==0;
k
		for kk=1:16
			if muus(kk)<0.001
				muus(kk)=0.001;
			end
		end
		xius(1+k/100)=xiauxus/100;
		xiauxus=0;
		mu10us(k/100)=fmuus(11);

		for kk=1:16
			if musgn(kk)<0.001
				musgn(kk)=0.001;
			end
		end
		xisgn(1+k/100)=xiauxsgn/100;
		xiauxsgn=0;
		mu10sgn(k/100)=fmusgn(11);

		xiopt(1+k/100)=xiauxopt/100;
		xiauxopt=0;
		mu10opt(k/100)=muopt(11);
	end
	s=[sign(randn); s(1:N-1)];
end

