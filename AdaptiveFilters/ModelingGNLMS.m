%
%	Modeling
%
%
% Last updated on April 28, 1998
%

% itn=input('\n No. of iterations?      ');
% sigman2=input('\n Variance of the plant noise?      ');
% sigman=sqrt(sigman2);
% wo=input('\n Plant impulse response (vector, w_o)?      ');
% a=size(wo);
% if a(1)<a(2)
%    wo=wo';
% end
% 
% N=input('\n Length of the model (N)?      ');
% 
% h=input('\n Coloring filter impulse response (vector, h)?     ');
% a=size(h);
% if a(1)<a(2)
%    h=h';
% end
% 
% Misad=input('\n Misadjustment (e.g., 0.1 for 10%) ?     ');
% mu=Misad/(N*(h'*h));
% 
% a=input('\n Do you wish to see the values of \n      eigenvalue spread, expected MSE, ... (Y/N)?      ','s');
% if (a=='y')|(a=='Y')
%    MMSE=sigman2;
%    R=corlnm2(h,N);
%    lambda=eig(R);
%    eignsprd=max(lambda)/min(lambda);
%    taumax=1/(4*mu*min(lambda));
%    taumin=1/(4*mu*max(lambda));
%    MSEaprx=MMSE*(1+Misad);
%    disp(' ')
%    disp(' ')
%    disp([' Eigenvalue spread = ' num2str(eignsprd)])
%    disp([' Maximum time constant of the learning curve = ' num2str(taumax)])
%    disp([' Minimum time constant of the learning curve = ' num2str(taumin)])
%    disp([' Expected steady-state MSE = ' num2str(MSEaprx)])
% end
% 
% runs=input('\n \n No. of runs (for ensemble averaging)? ');
% xi=zeros(itn,1);

itn=1000;
sigman2=0.0001; sigman=sqrt(sigman2);
N=20; M=5;
wo=randn(N,1);wo=wo/sqrt(wo'*wo);
h=0.9.^[0:20]'; h=h/sqrt(h'*h);
Misad=0.1;
mu=Misad/(N*(h'*h));
runs=10;
xi=zeros(itn,1);
xiGNLMS=zeros(itn,1);
mu=1; psi=1e-5;

for k=1:runs
	x=filter(h,1,randn(itn,1));
	d=filter(wo,1,x)+sigman*randn(itn,1);
	w=zeros(N,1);
    wGNLMS=w;
    X=zeros(N,M);
	for n=N+M:itn;
		xtdl=x(n:-1:n-N+1);
		e=d(n)-w'*xtdl;
		w=w+mu*e*xtdl/(xtdl'*xtdl);
		xi(n)=xi(n)+e^2;
        %%
        % GNLMS
        for k=1:M
            X(:,k)=x(n-k+1:-1:n-k+1-N+1);
        end
        e=d(n:-1:n-M+1)-X'*wGNLMS;
        %if rem(n,M)==0
            wGNLMS=wGNLMS+mu*X*inv(X'*X+psi*eye(M))*e; %end
        xiGNLMS(n)=xiGNLMS(n)+e(1)'*e(1);
        %%
	end
end
xi=xi/runs;
xiGNLMS=xiGNLMS/runs;

figure,semilogy(xi)
xlabel('NO. OF ITERATIONS')
ylabel('MSE')
hold on
semilogy(xiGNLMS,'r')
hold off
