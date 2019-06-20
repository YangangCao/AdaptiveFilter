%
%	Modeling (Standard RLS Algorithm I - Table 12.1)
%
%
% Last updated on April 28, 1998
%

itn=input('\n No. of iterations?      ');
sigman2=input('\n Variance of the plant noise?      ');
sigman=sqrt(sigman2);
wo=input('\n Plant impulse response (vector, w_o)?      ');
a=size(wo);
if a(1)<a(2)
   wo=wo';
end

N=input('\n Length of the model (N)?      ');

h=input('\n Coloring filter impulse response (vector, h)?     ');
a=size(h);
if a(1)<a(2)
   h=h';
end

lambda=input('\n Forgetting factor (lambda)?     ');
delta=input('\n Parameter delta?     ');

runs=input('\n \n No. of runs (for ensemble averaging)? ');
xi=zeros(itn,1);

for k=1:runs
	x=filter(h,1,randn(itn,1));
	d=filter(wo,1,x)+sigman*randn(itn,1);
	w=zeros(N,1);
   xtdl=zeros(size(w));
	Psi_inv=(1/delta)*eye(N);
	for n=1:itn
		xtdl=[x(n);xtdl(1:length(xtdl)-1)];
		u=Psi_inv*xtdl;
		k=u/(lambda+xtdl'*u);
		yhat=w'*xtdl;
		e=d(n)-yhat;
		w=w+k*e;
		Psi_inv=(1/lambda)*(Psi_inv-k*(xtdl'*Psi_inv));
		xi(n)=xi(n)+e^2;
	end
end
xi=xi/runs;
semilogy(xi)

