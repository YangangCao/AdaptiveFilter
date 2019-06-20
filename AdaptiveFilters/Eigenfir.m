%*******************************************************
%** This function can be added to MATLAB's vocabulary **
%** for designing the propotype filter of an M-band   **
%** complementary filter bank with the controlable    **
%** delay. The input parameters are:                  **
%**    M     ---  No. of bands of the filter bank;    **
%**    N     ---  Length of the designed filter;      **
%**    alpha ---  Roll-off factor;                    **
%**    del   ---  Group delay, i.e., Delta.           **
%** The output vector h stores the coefficients of    **
%** the designed FIR filter.                          **
%*******************************************************         
%
%   Format:       h = eigenfir(M,N,alpha,del)
% 
%
% Last updated on April 28, 1998
%

function h = eigenfir(M,N,alpha,del)
  
  ws = (1+alpha)/M;  
  corr_row = (-1)*ws*sinc(ws*[0:N-1]);
  corr_row(1) = 1-ws;
  corr = toeplitz(corr_row);

  left0 = [del-M:-M:1]; left0 = left0(size(left0,2):-1:1); 
  right0 = [del+M:M:N]; 
  zer0 = [left0 right0];
  corr(zer0,:) = []; corr(:,zer0) = [];

  [vect,eigen] = eig(corr); 
  [mi,mcnu] = min(diag(real(eigen)));
  vh = vect(:,mcnu); 
  j = 1; k = 1; 
  for i = 1:N
     if (i == zer0(j)) 
       h(i) = 0; 
       if j<size(zer0,2)
	 j = j+1;
       end 
     else
       h(i) = vh(k); 
       k = k+1; 
     end; 
  end; 
  h = h/h(del)/M;

