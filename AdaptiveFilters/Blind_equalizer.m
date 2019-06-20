%
%	Channel equalization
%
%
% Last updated on April 28, 1998
%

itn=2000;
ston=30;
N=15;
h=([0.3 1 -0.4]+1i*[-0.2 0 0.1])';
delay=8;
sigman=0;%10^(-ston/20)*sqrt(h'*h);
runs=1;
gamma=sqrt(2);

%% the case of l=1
xi=zeros(itn,1);
mu1=0.01;

for k=1:runs

	w=zeros(N,1);w(delay)=1;
	s=sign(randn(itn,1))+1i*sign(randn(itn,1));
	x=filter(h,1,s)+sigman*(randn(itn,1)+1i*randn(itn,1));

	for n=N:itn;
		xtdl=x(n:-1:n-N+1);
        y(n)=w'*xtdl;
		w=w-mu1*(conj(y(n))/abs(y(n)+0.0001))*(abs(y(n))-gamma)*xtdl;
		xi(n)=xi(n)+abs(abs(y(n))-gamma)^2;
	end
end
xi=xi/runs;

semilogy(xi)
xlabel('NO. OF ITERATIONS')
ylabel('MSE')
