%Lattice Joint Process Estimator.
%
%  This function program follows Table 11.5.
%
%  It is called by 'ltc_mdlg.m'.
%
% Last updated on April 28, 1998
%

function [kappa,c,b,e,P]=ljpe(kappa,c,x,d,b,P,mupo,muco,epsilon,beta);
N=length(c);
f=zeros(N,1);
bnew=f;
f(1)=x;
bnew(1)=x;
P(1)=beta*P(1)+0.5*(1-beta)*(f(1)^2+b(1)^2);
for m=2:N
	f(m)=f(m-1)-kappa(m-1)*b(m-1);
	bnew(m)=b(m-1)-kappa(m-1)*f(m-1);
	kappa(m-1)=kappa(m-1)+(2*mupo/(P(m-1)+epsilon))...
		*(f(m-1)*bnew(m)+b(m-1)*f(m));
	P(m)=beta*P(m)+0.5*(1-beta)*(f(m)^2+b(m)^2);
end
b=bnew;
y=c'*b;
e=d-y;
c=c+2*muco*e*(diag(P+epsilon)\b);


