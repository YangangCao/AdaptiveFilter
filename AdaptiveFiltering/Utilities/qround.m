function Vq = qround(V,b)
%QROUND Quantizes (rounds) decimal part.
%
%   Vq=QROUND(V,b) quantizes the decimal part of each element in
%                  vector V by rounding it to b bits.

	Vq=2^(-b)*round(V*2^b);
end