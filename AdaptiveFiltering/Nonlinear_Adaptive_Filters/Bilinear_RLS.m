%   Bilinear_RLS.m
%       Implements the Bilinear RLS algorithm for REAL valued data.
%       (Algorithm 11.3 - book: Adaptive Filtering: Algorithms and Practical
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
   x=Sx*(randn(dim,1));         % input signal
   xl(1)=0; xl(2:dim)=x(1:dim-1);
   w=zeros(4,dim);           % initial coefficient vector
   d=zeros(dim,1);
   yl(1)=0; dl(1)=0;
   
   Sd=eye(4);
   for i=1:dim
      d(i)=-.6*dl(i)+x(i)+.01*x(i)*dl(i)+.02*xl(i)*dl(i)+n(i); % unknown system output sample
      uxl(:,i)=[x(i) dl(i) x(i)*dl(i) xl(i)*dl(i)]'; % new input vector
      elinha(i)=d(i)-w(:,i)'*uxl(:,i)+n(i);          % error sample
      psi=Sd*uxl(:,i);
      Sd=(1/lambda)*(Sd-(psi*psi')/(lambda+psi'*uxl(:,i)));
      w(:,i+1)=w(:,i)+elinha(i)*Sd*uxl(:,i);         % new coefficient vector
      y(i)=w(:,i+1)'*uxl(:,i);                       % output sample
      yl(i+1)=yl(i);
      dl(i+1)=d(i);
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


