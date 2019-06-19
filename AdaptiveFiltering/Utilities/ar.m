function a1=ar(n,Cd)
%AR obtains a first-order AR process with given characteristics.
%
%   a1=AR(n,Cd) determines the coefficient a1 of a first-order AR
%               process such that, when this is applied at the 
%               input of an order-n tapped delay line, the 
%               eigenvalue spread of the correlation matrix of 
%               the input vector is Cd.

   min=-1;	% minimum coefficient
   max=0;	% maximum coefficient
   C=1;		% initial (minimum) eigenvalue spread

   while abs((C-Cd)/Cd)>1e-6,	% while difference < tolerance
      a1=(min+max)/2;		% mean coefficient
      for i=1:n+1,
         for j=1:n+1,
            A(i,j)=(-a1)^abs(i-j);	% new correlation matrix
         end
      end
      C=cond(A);		% new eigenvalue spread
      if C<Cd,			% if spread low, try a1 lower
         max=a1;
      else 			% if spread high, try a1 higher
         min=a1;
      end
   end
end