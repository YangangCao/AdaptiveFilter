
%LMS2 Problem 1.1.1.1.2.1
%
%   'ifile.mat' - input file containing:
%      I - members of ensemble
%      K - iterations
%      sigmax - standard deviation of input
%      Wo - coefficient vector of plant
%      sigman - standard deviation of measurement noise
%      mu - convergence factor
%      b - bits in decimal part
% 
%   'ofile.mat' - output file containing:
%      ind - sample indexes 
%      MSE - mean-square error
%      MSNDW - mean-square norm of coefficient-error vector

clear all		% clear memory
load ifile;		% read input variables
L=length(Wo);		% plant and filter length
N=L-1;			% plant and filter order
MSE=zeros(K,1);		% prepare to accumulate MSE*I
MSNDW=zeros(K,1);	% prepare to accumulate MSNDW*I

for i=1:I,		% ensemble
   X=zeros(L,1);    	% initial memory
   W=zeros(L,1);	% initial coefficient vector
   x=randn(K,1)*sigmax;		% input 
   n=randn(K,1)*sigman;		% measurement noise 
   for k=1:K,		% iterations
      X=[x(k)		
         X(1:N)];	% new input vector
      d=Wo'*X;		% desired signal sample
      y=W'*X;		% output sample
      e=d+n(k)-y;	
      e=qround(e,b);	% error sample
      W=W+2*mu*e*X;	
      W=qround(W,b); 	% new coefficient vector
      MSE(k)=MSE(k)+e^2;	% accumulate MSE*I
      MSNDW(k)=MSNDW(k)+norm((Wo-W),2)^2; 
                        % accumulate MSNDW*I
   end
end

ind=0:(K-1);		% sample indexes
MSE=MSE/I;		% calculate MSE
MSNDW=MSNDW/I;		% calculate MSNDW
save ofile ind MSE MSNDW;	% write output variables
