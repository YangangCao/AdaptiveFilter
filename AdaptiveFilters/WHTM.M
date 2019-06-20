function T=whtm(n)
%This function produces a matrix of dimension n by n
%and with elements of the WHT.
%
%Format:			T=whtm(n)
%
%NOTE: n must be a power of 2.
%
p=round(log(n)/log(2));
T=ones(n)/sqrt(n);
for m=0:n-1
	for k=0:n-1
		count=0;
		for l=0:p-1
			count=count+(fix(m/(2^l))-2*fix(m/(2^(l+1)))) ...
					*(fix(k/(2^(p-1-l)))-2*fix(k/(2^(p-l))));
		end
		T(m+1,k+1)=T(m+1,k+1)*(-1)^count;
	end
end
