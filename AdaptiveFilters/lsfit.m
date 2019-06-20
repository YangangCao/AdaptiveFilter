%
%This function finds the coefficients of an ARMA model which
%for a given input sequence x(n) results in an output which
%matches best the sequence d(n) in the least-squares sence.
%
%  This program is called by 'iirdsgn.m'.
%
%
% Last updated on April 28, 1998
%

function [A,B]=lsfit(x,d,N_d,N_n);
N=length(x);
C=zeros(N-N_n-N_d+1,N_n+N_d);
for k=N_n+N_d:N
	C(k-N_n-N_d+1,:)=[x(k:-1:k-N_n+1)' d(k-1:-1:k-N_d)'];
end
BA=(C'*C)\(C'*d(N_n+N_d:N));
B=BA(1:N_n);
A=[1; -BA(N_n+1:N_n+N_d)];



