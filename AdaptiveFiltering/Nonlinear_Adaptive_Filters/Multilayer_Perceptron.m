%   Multilayer_Perceptron.m
%       Implements the Multilayer Perceptron algorithm for REAL valued data.
%       (Algorithm 11.4 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, Diniz)
%  
%   Input parameters:
%       Nr     : members of ensemble.
%       dim    : iterations.
%       Sx     : standard deviation of input.
%       Sn     : standard deviation of measurement noise.
%       lambda : exponential weighting factor.
%       Nneur  : number of neurons
%       u      : convergence factor
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
u      = 0.01;


% Body:
for j=1:Nr
   n=Sn*randn(dim,1);           % noise at system output   
   x=Sx*(randn(dim,1));            % input signal 
   xl(1)=0; xl(2:dim)=x(1:dim-1);
   d=zeros(dim,1);
   yl(1)=0; dl(1)=0;
   % initial training vectors
   w1=.2*randn(Nneur,3); w2=.2*randn(Nneur,Nneur); w3=.2*randn(Nneur,1);
   b1=.1*randn(Nneur,1); b2=.1*randn(Nneur,1);  b3=.1*randn(1,1);
   for i=1:dim
      d(i)=-.6*dl(i)+x(i)+.1*x(i)*dl(i)+.2*xl(i)*dl(i)+n(i); % unknown system output sample
      uxl=[x(i) dl(i) xl(i)]';     % new input vector
      aux=w2*sgm(w1*uxl-b1,1,2)-b2;
      dr(i)=sgm(aux,1,2)'*w3-b3;   % neural network output sample
      er_out=d(i)-dr(i);           % error sample
      er_hid2=er_out*w3.*sgd(aux,1,2);
      er_hid=w2'*er_hid2.*sgd(w1*uxl-b1,1,2);
     dw3=2*u*er_out*sgm(aux,1,2);
     db3=-2*u*er_out;
     dw2=2*u*er_hid2*sgm(w1*uxl-b1,1,2)';
     db2=-2*u*er_hid2;
     dw1=2*u*er_hid*uxl';
     db1=-2*u*er_hid;       
      w1=w1+dw1; b1=b1+db1;        % new coefficient vectors for layer 1
      w2=w2+dw2; b2=b2+db2;        % new coefficient vectors for layer 2
      w3=w3+dw3; b3=b3+db3;        % new coefficient vectors for layer 3
      mse(j,i)=(d(i)-dr(i))^2;
      dl(i+1)=d(i);
   end
end  

MSE=mean(mse);


% Output:
figure,
plot(10*log10(MSE));
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

