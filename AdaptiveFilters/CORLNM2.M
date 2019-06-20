function R=corlnm2(h,N)
%This function obtains the correlation matrix of the output
%of an FIR filter which is derived by a white noise, and its
%coefficients are given in vector "h".
%"N" specifies the dimensions of "R".
%
%	R=corlnm2(h,N)
j=sqrt(-1);
a=size(h);
if a(1)>a(2)
   h=h';
end
m=length(h);

for k=1:N
	if k<=m
		R(k)=h(k:m)*h(1:m-k+1)';
	else	R(k)=0; end
end
R=toeplitz(R);
