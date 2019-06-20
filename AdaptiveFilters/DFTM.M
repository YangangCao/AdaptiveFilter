function T=dftm(n)
%This function produces the coeffcients of a DFT matrix
%of dimension n by n.
%
%Format:   T=dftm(n)
%
T=zeros(n);
i=sqrt(-1);
a=exp(i*2*pi/n);
for m=1:n
	for k=1:n
		T(m,k)=a.^((m-1).*(k-1))./sqrt(n);
	end
end
