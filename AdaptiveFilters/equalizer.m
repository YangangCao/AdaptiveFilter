%
%	Channel equalization
%
%
% Last updated on April 28, 1998
%

itn=input('\n No. of iterations?      ');
ston=input('\n Signal-to-noise ratio (dB)?      ');
N=input('\n Length of equalizer (N)?      ');
h=input('\n Channel response (vector, h)?     ');
a=size(h);
if a(1)<a(2)
   h=h';
end

Misad=input('\n Misadjustment (e.g., 0.1 for 10%) ?     ');
delay=input('\n Delay (Delta)?     ');
sigman=10^(-ston/20)*sqrt(h'*h);
mu=Misad/(N*(h'*h)*(1+10^(-ston/10)));

a=input('\n Do you wish to see the values of \n      eigenvalue spread, expected MSE, ... (Y/N)?      ','s');
if (a=='y')|(a=='Y')
   R=corlnm2(h',N)+sigman^2*eye(N);
   p=[zeros(delay-length(h)+1,1); h(length(h):-1:1); zeros(N-delay-1,1)];
   MMSE=1-p'*inv(R)*p;
   lambda=eig(R);
   eignsprd=max(lambda)/min(lambda);
   taumax=1/(4*mu*min(lambda));
   taumin=1/(4*mu*max(lambda));
   MSEaprx=MMSE*(1+Misad);
   disp(' ')
   disp(' ')
   disp([' Eigenvalue spread = ' num2str(eignsprd)])
   disp([' Maximum time constant of the learning curve = ' num2str(taumax)])
   disp([' Minimum time constant of the learning curve = ' num2str(taumin)])
   disp([' Expected steady-state MSE = ' num2str(MSEaprx)])
end

runs=input('\n \n No. of runs (for ensemble averaging)? ');

xi=zeros(itn,1);
mu=Misad/(N*(h'*h)*(1+10^(-ston/10)));

for k=1:runs

	w=zeros(N,1);
	s=sign(randn(itn,1));
	x=filter(h,1,s)+sigman*randn(itn,1);

	for n=N:itn;
		xtdl=x(n:-1:n-N+1);
		e=s(n-delay)-w'*xtdl;
		w=w+2*mu*e*xtdl;
		xi(n)=xi(n)+e^2;
	end
end
xi=xi/runs;

semilogy(xi)
xlabel('NO. OF ITERATIONS')
ylabel('MSE')
