function T=dhtm(n)
%This function produces a matrix of dimension n by n
%and with elements of the discrete Hartley transform (DHT).
%
%Format:   T=dhtm(n)
%
T=zeros(n);
for m=1:n
	for k=1:n
		teta=2*pi*(m-1)*(k-1)/n;
		T(m,k)=sqrt(1/n)*(cos(teta)+sin(teta));
	end
end
