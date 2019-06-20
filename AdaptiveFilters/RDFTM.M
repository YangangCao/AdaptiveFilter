function T=rdftm(n)
%This function produces the coeffcients of a RDFT matrix
%of dimension n by n.
%
%Format:   T=rdftm(n)
%
T=zeros(n);
T(1,:)=ones(1,n)/sqrt(n);
T(2:n/2,:)=sqrt(2/n)*cos((2*pi/n)*[1:n/2-1]'*[0:n-1]);
T(n/2+1,:)=(1/sqrt(n))*(-1).^[0:n-1];
T(n/2+2:n,:)=sqrt(2/n)*sin((2*pi/n)*[n/2+1:n-1]'*[0:n-1]);
