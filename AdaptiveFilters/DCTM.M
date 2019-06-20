function T=dctm(n)
%This function produces a matrix of dimension n by n
%and with elements of the DCT.
%
%Format:   T=dctm(n)
%
T=zeros(n);
for m=1:n
	for k=1:n
		if m==1
			T(m,k)=1./sqrt(n);
		else
			T(m,k)=sqrt(2./n).*cos(pi.*(2.*k-1) ...
						 .*(m-1)./(2.*n));
								end
	end
end

