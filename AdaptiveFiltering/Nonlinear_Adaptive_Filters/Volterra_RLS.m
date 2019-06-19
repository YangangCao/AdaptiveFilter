%   Volterra_RLS.m
%       Implements the Volterra RLS algorithm for REAL valued data.
%       (Algorithm 11.2 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, Diniz)
%  
%   Input parameters:
%       Nr     : members of ensemble.
%       dim    : iterations.
%       Sx     : standard deviation of input.
%       Sn     : standard deviation of measurement noise.
%       lambda : exponential weighting factor.
%      
%   Output parameters:
%       MSE    : mean-square error.
%  
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com & guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom           & markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com     & wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com          & wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                             diniz@lps.ufrj.br
%


clear all;		% clear memory


% Input: 
Nr     = 100;
dim    = 5e2;
Sx     = 1; 
Sn     = 1e-1;
lambda = 0.98;


% Body:
for j=1:Nr
   n=Sn*randn(dim,1);        % noise at system output 
   x=Sx*randn(dim,1);       % input signal 
   xl1=zeros(dim,1); xl2=xl1;
   xl1(2:dim)=x(1:dim-1); xl2(3:dim)=x(1:dim-2);
   d=zeros(dim,1);
   d=-.76*x-xl1+xl2+.5*x.^2+2*x.*xl2-1.6*xl1.^2+1.2*xl2.^2+.8*xl1.*xl2+n; ...
     % unknown system output
   w=zeros(9,dim);           % initial coefficient vector
   uxl=[x xl1 xl2 x.^2 x.*xl1 x.*xl2 xl1.^2 xl1.*xl2 xl2.^2]'; % input vectors
   
   Sd=eye(9);
   for i=1:dim
      elinha(i)=d(i)-w(:,i)'*uxl(:,i)+n(i);  % error sample
      psi=Sd*uxl(:,i);
      Sd=(1/lambda)*(Sd-(psi*psi')/(lambda+psi'*uxl(:,i)));
      w(:,i+1)=w(:,i)+elinha(i)*Sd*uxl(:,i); % new coefficient vector
      y(i)=w(:,i+1)'*uxl(:,i);               % output sample
      e(i)=d(i)-y(i)+n(i);
   end
   mse(j,:)=e.^2;
end

MSE=mean(mse);


% Output:
figure,
plot(10*log10(MSE));
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

