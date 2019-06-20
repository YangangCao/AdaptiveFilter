%
%	Antenna Array (Real-valued signal at an IF are considered)
%
%
% Last updated on April 28, 1998
%

itn=input('\n No. of iterations?      ');
omegao=input('\n angular frequency omega_o (rad/s)?     ');
thetao=input('\n angle of arrival, theta_o (degrees)?     ');
thetao=pi*thetao/180;
deltao=pi*sin(thetao);
Misad=input('\n Misadjustment (e.g., 0.1 for 10%) ?     ');
sigma_a=input('\n Variance of the desired signal, sigma_alpha?     ');
sigma_b=input('\n Variance of jammer signal, sigma_beta?     ');
traceR=2*(0.5*sigma_a+0.5*sigma_b);
mu=Misad/traceR;
a=input('\n Do you wish to see the values of \n      eigenvalue spread, expected MSE, ... (Y/N)?      ','s');
if (a=='y')|(a=='Y')
   R=diag(traceR*[0.5 0.5]);
   p=0.5*[sigma_a+sigma_b*cos(deltao); sigma_b*sin(deltao)];
   MMSE=0.5*traceR-p'*(R\p);
   MSEaprx=MMSE*(1+Misad);
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

xi=zeros(itn,1);
runs=input('\n \n No. of runs (for ensemble averaging)? ');

for k=1:runs
	a=sqrt(sigma_a)*randn(itn,1);	%desired signal
	b=sqrt(sigma_b)*randn(itn,1);	%jammer
   theta_a=2*pi*rand;	         %randon phase for desired signal
	theta_b=2*pi*rand;	         %random phase for jammer
	%
	% "xp" is the signal picked up at the primary.
	% 
	xp=a.*cos([1:itn]'*omegao+theta_a)+b.*cos([1:itn]'*omegao-deltao+theta_b);
	%
	% "x1" and "x2" are the signals at the reference taps.
	%
	x1=a.*cos([1:itn]'*omegao+theta_a)+b.*cos([1:itn]'*omegao+theta_b);
	x2=a.*sin([1:itn]'*omegao+theta_a)+b.*sin([1:itn]'*omegao+theta_b);
	w=[0 0]';

	for n=1:itn
		xtdl=[x1(n); x2(n)];
		e=xp(n)-w'*xtdl;
		w=w+2*mu*e*xtdl;
		xi(n)=xi(n)+e^2;
	end
end
xi=xi/runs;
figure(1); semilogy(xi)
theta=2*pi*[0:0.01:1];
gain=(cos(pi*sin(theta))-w(1)).^2+(sin(pi*sin(theta))-w(2)).^2;
figure(2);
f=polar(theta,gain);
set(f,'LineWidth',1.5)
