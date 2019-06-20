function h = lorentz(N,L,t50,delay)
% 	h = lorentz(N,L,t50,delay)
%
% To compute the Lorentzian pulse shape.
%
%  N = total number of samples (odd number)
%  L = over-sampling factor
%  D = t50 = half height width (T = 1sec.) 
%  delay = shift of the pulse center to right.
%
% Last updated on April 28, 1998
%
  t=2*([-(N-1)/2:(N-1)/2]+delay*L)/(L*t50);
  h = (1 + t.^2).^(-1); 


