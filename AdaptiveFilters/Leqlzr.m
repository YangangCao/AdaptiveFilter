function w=leqlzr(h,L,delta,N,d,sigman2);
%
%This functin calculate the coefficient of a linear fractionally
%tap-spaced equalizer.
%
%	It has the following format:
%
%		w=leqlzr(h,L,delta,N,d,sigman2);
%	where
%		"h" is the channel response at the rate of T/L
%		"T" is the bit interval
%		"L" is the oversampling ratio
%		"delta" is the equalizer response delay in the units of T/L
%			before its first non-zero sample.
%		"N" is the equalizer length
%		"d" is the desired response (starting with a non-zero sample).
%		"sigmav2" is the noise variance.
% 
%  This function is called by 'iirdsgn.m'.
%
% Last updated on April 28, 1998
%

M=length(h);
P=ceil((N+M)/L);
H=zeros(P,N);
h=[zeros(1,N+M-1-rem(delta,L)) h zeros(1,N+M-1)];
d=[zeros(1,floor(delta/L)) d];
d=[d zeros(1,P-length(d))]';
for k=1:P
	H(k,:)=h(N+M+(k-1)*L:-1:(N+M+(k-1)*L-(N-1)));
end
w=(H'*H+sigman2*eye(N))\(H'*d);
