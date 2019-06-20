function T=dstm(n)
%This function produces a matrix of dimension n by n
%and with elements of the DST of type I.
%
%Format:   T=dstm(n)
%
T=zeros(n);
for m=1:n
	for k=1:n
		T(m,k)=sqrt(2/(n+1))*sin(m*k*pi/(n+1));
	end
end
