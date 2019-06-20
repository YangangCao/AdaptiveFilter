%Fast LMS-Newton algorithm when applied to the modeling of an
%unknown plant: Algorithm 2, Table 11.7.
%
%  
%  Lattice structures are used to implement the forward 
%  and backward filters in the AR modeling part of the algorithm.
%
%  The common input to the plant and adaptive filter is generated
%  by passing a white noise sequence through a ploe-zero (ARMA)
%  transfer function.
%
itn=input('\n No. of iterations?      ');
N=input('\n Length of the plant/adaptive filter (N)?      ');
M=input('\n Order of AR model (M)?     ');
B=input('\n Noise coloring filter (numerator)?     ');
A=input('\n Noise coloring filter (denominator)?     ');
dummy=input('\n Do you wish to select the plant rondomly (Y/N)?      ','s');
if (dummy=='y')|(dummy=='Y')
   wo=randn(N,1);
   wo=wo/sqrt(wo'*wo);
else
   wo=input('\n Plant impulse response (vector, w_o)?      ');
   dmy=size(wo);
   if dmy(1)<dmy(2)
      wo=wo';
   end
end
sigman2=input('\n Variance of the plant noise?      ');
sigman=sqrt(sigman2);
Misad=input('\n Misadjustment (e.g., 0.1 for 10%) ?     ');
mu=Misad/N;
mu_po=input('\n Step-size parameter mu_po?     ');
alpha=input('\n Parameter alpha?     ');
beta=input('\n Parameter beta?     ');
epsilon=input('\n Parameter epsilon?     ');
runs=input('\n \n No. of runs (for ensemble averaging)? ');
xi=zeros(itn,1);

for k=1:runs
   xin=randn(itn,1);
   xin=filter(B,A,xin);
   d=filter(wo,1,xin)+sigman*randn(itn,1);  %Plant output
   P=ones(M+1,1);
   w=zeros(N,1);
   kappa=zeros(M,1);
   b=xin(N+M:-1:N);
   bp=zeros(M,1);           %'bp' denotes the vector b' of Table 11.7
                            %Similarly 'fp' is used below for f'. 
   u_a=zeros(N,1);

   for n=N+M+1:itn
   %
	%	Lattice Predictor
	%
   b_old=b;
   f(1)=xin(n);
   b(1)=f(1);
   P(1)=beta*P(1)+0.5*(1-beta)*(f(1)^2+b_old(1)^2);
   for m=1:M
      f(m+1)=f(m)-kappa(m)*b_old(m);
      b(m+1)=b_old(m)-kappa(m)*f(m);
      kappa(m)=kappa(m)+2*mu_po*(P(m)+epsilon)^(-1)*(f(m)*b(m+1)+b_old(m)*f(m+1));
      P(m+1)=beta*P(m+1)+0.5*(1-beta)*(f(m+1)^2+b_old(m+1)^2);
   	if abs(kappa(m))>alpha
			kappa(m)=sign(kappa(m))*alpha;
		end
	end
   %
   %  u_a(n) update
   %
   u_a(2:N)=u_a(1:N-1);
   b_old=bp;
   fp(1)=b(M+1);
   bp(1)=fp(1);
   for m=1:M
      fp(m+1)=fp(m)-kappa(m)*b_old(m);
      bp(m+1)=b_old(m)-kappa(m)*fp(m);
	end
   u_a(1)=fp(M+1)/(P(M+1)+epsilon);
   %
   %   Filtering
   %
   y=w'*xin(n-M:-1:n-N-M+1);
   e=d(n-M)-y;
	w=w+2*mu*e*u_a;
	xi(n)=xi(n)+e*e;
   end
end
xi=xi/runs;
n=[1:itn];
semilogy(n,xi)


