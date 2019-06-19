%   Complex_Radial_Basis_Function.m
%       Implements the Complex Radial Basis Function algorithm for COMPLEX valued data.
%       (Algorithm 11.6 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, Diniz)
%  
%   Input parameters:
%       Nr     : members of ensemble.
%       dim    : iterations.
%       Sx     : standard deviation of input.
%       Sn     : standard deviation of measurement noise.
%       Nneur  : number of neurons.
%       ur     : convergence factor for reference vector.
%       uw     : convergence factor for coefficient vector.
%       us     : convergence factor for spread.
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
Nneur  = 10;
ur     = 0.01;
uw     = 0.01;
us     = 0.01;


% Body:
for j=1:Nr
   n=Sn*randn(dim,1);           % noise at channel output 
   x=Sx*randn(dim,1)+sqrt(1)*Sx*randn(dim,1);     % input signal
   xl1=zeros(dim,1); xl2=xl1;
   xl1(1)=0; xl1(2:dim)=x(1:dim-1);
   xl2(1)=0; xl2(2:dim)=xl1(1:dim-1); 
   d=-.08*x-.15*xl1+.14*xl2+.055*x.^2+.3*x.*xl2-.16*xl1.^2+.14*xl2.^2+n; ...
     % unknown system output
   w=randn(Nneur,dim); w(:,1)=randn(Nneur,1);  % initial coefficient vector
   vet=.5*randn(Nneur,3);                      % initial reference vector
   sigma=ones(Nneur,1);                        % initial neuron spread
   for i=1:dim
      uxl(:,i)=[x(i) xl1(i) xl2(i)].';      % new input vector
      dis=dist(uxl(:,i).',vet.').';
      fdis=exp(-(dis.^2)./(sigma.^2));
      e(i)=d(i)-w(:,i)'*fdis;              % error sample
      w(:,i+1)=w(:,i)+2*uw*e(i)*fdis;      % new coefficient vector
      sigma=sigma+2*us*fdis.*(real(e(i)).*real(w(:,i))+imag(e(i)).*imag(w(:,i))).*(dis.^2)./(sigma.^3); % new spread 
      for p=1:Nneur
        vet(p,:)=vet(p,:)+2*ur*fdis(p)*(real(e(i))*real(w(p,i))*real(uxl(:,i).'-vet(p,:))+sqrt(-1)*imag(e(i))*imag(w(p,i))*imag(uxl(:,i).'-vet(p,:)))/ ...
		 (sigma(p)^2); % new reference vector
     end
     y(i)=w(:,i)'*fdis;                    % output sample
   end
   mse(j,:)=e.^2;
end

MSE=mean(mse);


% Output:
figure,
plot(10*log10(MSE));
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');


