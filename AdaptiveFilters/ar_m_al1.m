%Fast LMS-Newton algorithm when applied to the modeling of an
%unknown plant: Algorithm 1, Table 11.6.
%
%  
%  Lattice structures are used to implement the backward filter 
%  in the AR modeling part of the algorithm.
%  The Levinson algorithm of Table 11.1 is used to convert the PARCOR
%  coefficients of the lattice to the corresponding transversal
%  structure coefficients (a_{i,j}'s).
%
%  The common input to the plant and adaptive filter is generated
%  by passing a white noise sequence through a ploe-zero (ARMA)
%  transfer function.
%
itn=input('\n No. of iterations?      ');
N=input('\n Length of the plant/adaptive filter (N)?      ');
M=input('\n Order of AR model (M)?     ');
B_c=input('\n Noise coloring filter (numerator)?     ');
A_c=input('\n Noise coloring filter (denominator)?     ');
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
   xin=filter(B_c,A_c,xin);
   d=filter(wo,1,xin)+sigman*randn(itn,1);  %Plant output
   P=ones(M+1,1);
   Pbar=ones(N+M,1);
   w=zeros(N,1);
   kappa=zeros(M,1);
   f=ones(M,1)*xin(N+M);
   b=xin(N+M:-1:N);
   bb=zeros(N,1);
   u=bb;
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
   bb=[b; bb(M+1:N-1)];     %un-normalized backward errors 
                  %normalization of these to get b_h and b_t is done later.
	%
	%	Conversion of a_M,i coefficients (Levinson Algorithm).
   %
   %  This part of the program builds up the vector Pbar
   %  and the matrix
   %
   %         |a_{1,1}    0       0    0 .......... 0|
   %         |a_{2,1} a_{2,2}    0    0 .......... 0|
   %         |a_{3,1} a_{3,2} a_{3,3} 0 .......... 0|
   %     A = |......................................|
   %         |......................................| 
   %         |a_{M,1} a_{M,2} a_{M,3} ...... a_{M,M}|
   %         
   %
   Pbar(1)=P(1);
   A=zeros(M);
   A(1,1)=kappa(1);
	for m=1:M-1
   	Pbar(m+1)=Pbar(m)*(1-kappa(m)^2);
		for j=m:-1:1
			A(m+1,j)=A(m,j)-kappa(m+1)*A(m,m+1-j);
		end
      A(m+1,m+1)=kappa(m+1);          
   end
   Pbar(2*M+2:N+M)=Pbar(2*M+1:N+M-1); %Extended Pbar to be used for normalization 
   Pbar(M+1:2*M+1)=Pbar(M)*(1-kappa(M)^2)*ones(M+1,1);  %of b_h (see below).
                                         
   %
   %  Build L_tl 
   %
   L_tl=zeros(2*M+1,M+1);
   L_tl(1:M+1,:)=eye(M+1);
   for m=2:M+1
      L_tl(m,1:m-1)=-A(m-1,m-1:-1:1);
   end
   for m=M+2:2*M+1
      L_tl(m,m-M:M+1)=-A(M,M:-1:m-M-1);
   end
   %
   %  Build L_br
   %
   L_br=eye(M);
   for m=2:M
      L_br(m,1:m-1)=-A(M,m-1:-1:1);
   end
   %
   %  u(n) update
   %
   u(M+2:N-M)=u(M+1:N-M-1);
   b_h=bb(1:2*M+1)./Pbar(1:2*M+1);
   u(1:M+1)=L_tl'*b_h;
   b_t=Pbar(N+M)\bb(N-M+1:N);
   u(N-M+1:N)=L_br'*b_t;
   %
   %   Filtering
   %
   y=w'*xin(n:-1:n-N+1);
   e=d(n)-y;
	w=w+2*mu*e*u;
	xi(n)=xi(n)+e*e;
   end
end
xi=xi/runs;
n=[1:itn];
semilogy(n,xi)

