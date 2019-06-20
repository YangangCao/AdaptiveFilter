%
%   Modeling using lattice and transversal structures.
%   Comparison based on LMS algorithm.
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

B=input('\n Coloring filter (numerator, B)?     ');
A=input('\n Coloring filter (denominator, A)?     ');

mupo=input('\n PARCOR coefficients step-size parameter, mu_p,o?     ');
muco=input('\n Joint process estimator coefficients step-size parameter, mu_c,o?      ');
epsilon=input('\n Parameter epsilon?     ');
beta=input('\n Parameter beta?     ');

disp(' ')
disp(' Do you wish to stop adaptation of PARCOR Coefficients')
a=input(' after certain no. of iterations (Y/N)?      ','s');
s_PARCORs=itn+1;
if (a=='y')|(a=='Y')
   s_PARCORs=input('\n Iteration no. at which PARCORs adaption shall be stoped?     ');
end

xiltc=zeros(itn,1);
xilms=xiltc;
mu=muco;
runs=input('\n No. of runs (for ensemble averaging)? ');

for k=1:runs
   %
   %   Initialization
   %
   x=randn(itn,1);
   x=filter(B,A,x);
   x=x/std(x);              %Set the variance of input equal to one.
   d=filter(wo,1,x)+sigman*randn(size(x));
   kappa=zeros(N-1,1);
   c=zeros(N,1);
   b=zeros(N,1);
   P=ones(N,1);
   xin=zeros(N,1);
   w=zeros(N,1);
   %
   %      Iteration loop
   %
   for n=1:itn
      if n==s_PARCORs 
         mupo=0; 
      end
      %
      %   Lattice LMS update
      %
      [kappa,c,b,e,P]=ljpe(kappa,c,x(n),d(n),b,P,mupo,muco,epsilon,beta);
      xiltc(n)=xiltc(n)+e*e;
      %
      %   Transversal LMS update
      %
      xin=[x(n); xin(1:N-1)];
      y=w'*xin;
      e=d(n)-y;
      w=w+2*mu*e*xin;
      xilms(n)=xilms(n)+e*e;
   end
end

nn=1:itn;
semilogy(nn,xiltc,'-',nn,xilms,'--')
axis([1 itn 10^(-5) 10])
xlabel('NO. OF ITERATIONS')
ylabel('MSE')
title('Full-line: LMS-Lattice,   Dashed-line: Cinventional LMS')

