%
%	Transversal line enhancer
%
%  The line enhancer input consists of a sine-wave plus a white noise:
%
%                       x(n) = sin(omega_o*n)+v(n)
%
%  The amplitude of the sine-wave is set equal to one.
%  The variance of the noise, v(n), is thus determined
%  according to the specified value of the signal-to-noise ratio.
%
%
% Last updated on April 28, 1998
%

itn=input('\n No. of iterations?      ');
ston=input('\n Signal-to-noise ratio (dB)?      ');
sigman=(10^(-ston/20))/sqrt(2);	
N=input('\n Length of the line enhancer (N)?      ');
omega=input('\n the sine-wave frequency, omega_o?     ');
Misad=input('\n Misadjustment (e.g., 0.1 for 10%) ?     ');
R=toeplitz(0.5*cos(omega*[0:N-1]))+(sigman^2)*eye(N);
mu=Misad/trace(R);


a=input('\n Do you wish to see the values of \n      eigenvalue spread, expected MSE, ... (Y/N)?      ','s');
if (a=='y')|(a=='Y')
   p=0.5*cos(omega*[1:N]');
   MMSE=0.5+sigman^2-p'*inv(R)*p;
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

wflag=' ';
while (wflag~='r')&(wflag~='R')&(wflag~='z')&(wflag~='Z')
   disp(' ')
   disp(' ')
   disp(' How do you want the line enhancer tap weight to be initialized?')
   disp('         enter "z" to initialize to zero,')
   wflag=input('         enter "r" to initialize to some randon values:    ','s');
end

runs=input('\n \n No. of runs (for ensemble averaging)? ');
y=zeros(itn,1);
xi=zeros(itn,1);

for k=1:runs
	x=sin(omega*[1:itn]'+2*pi*randn)+sigman*randn(itn,1);
   if (wflag=='z')|(wflag=='Z')
      w=zeros(N,1);
   elseif (wflag=='r')|(wflag=='R')
      w=randn(N,1);
   end
   
   for n=N+1:itn;
		xtdl=x(n-1:-1:n-N);
		y(n)=w'*xtdl;
		e=x(n)-y(n);
		w=w+2*mu*e*xtdl;
		xi(n-N)=xi(n-N)+e^2;
	end
end
xi=xi/runs;
semilogy(xi)
xlabel('NO. OF ITERATIONS')
ylabel('MSE')


