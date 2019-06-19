%   olsblms.m      
%       Implements the Open-Loop Subband (LMS) Adaptive-Filtering Algorithm for REAL valued data.
%       (Algorithm 12.1 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, Diniz)
%  
%   Input parameters:
%       Nr     : members of ensemble.
%       N      : iterations.
%       Sx     : standard deviation of input.
%       Sn     : standard deviation of measurement noise.
%       B      : coefficient vector of plant (numerator).
%       A      : coefficient vector of plant (denominator).
%       u      : convergence factor.
%       gamma  : small constant to prevent the updating factor from getting too large.
%       a      : small factor for the updating equation for the signal energy in the subbands.
%       M      : number of subbands.
%       hk     : analysis filterbank, with the filter coefficients in the rows.
%       fk     : synthesis filterbank, with the filter coefficients in the rows.
%       Nw     : order of the adaptive filter in the subbands.
%      
%   Output parameters:
%       MSE    : mean-square error.
%  
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com & guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom           & markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com     & wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com          & wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                             diniz@lps.ufrj.br
%


clear all;		   % clear memory


% Input: 
load cosmod_4_64;          % load filter bank parameters: M, hk, fk
Nr     = 100;
N      = 5e2;
Sx     = 1; 
Sn     = 1e-1;
B      = [0.32,-0.3,0.5,0.2];
A      = 1; 
u      = 0.1;
gamma  = 1e-2;
L      = M;                % interpolation/decimation factor
a      = 0.01;
Nw     = 5;


for l=1:Nr
    disp(['Run: ' num2str(l)])
    nc=sqrt(Sn)*randn(1,N*M);     % noise at system output 
    xin=(rand(1,N*M)-.5);
    xin=sqrt(Sx)*sqrt(12)*xin;    % input signal
    d=filter(B,A,xin);            % desired signal
    d=d+nc;                       % unknown system output

    % Analysis of desired and input signals
    for k=1:M
	xaux=filter(hk(k,:),1,d);
	dsb(k,:)=xaux(find(mod((1:length(xaux))-1,L)==0)); % desired signal split in subbands
	xaux=filter(hk(k,:),1,xin);
	xsb(k,:)=xaux(find(mod((1:length(xaux))-1,L)==0)); % input signal split in subbands
    end;
  
    for m=1:M
	w_ol(m,:)=zeros(1,Nw+1);      % initial coefficient vector for
                                      % subband m
	x_ol(m,:)=[zeros(1,Nw+1)];    % initial input vector for subband m
	sig_ol(m)=0;
    end

    for k=1:N
      if mod(k,200)==0 disp(['Iteration: ' num2str(k)]), end;
      for m=1:M 
	  % Open-loop
	  % --------------------------------------------------
	  x_ol(m,:)=[xsb(m,k) x_ol(m,1:Nw)]; % new input vector for subband m
	  xaux=shiftdim(x_ol(m,:),1);
	  waux=shiftdim(w_ol(m,:),1);
	  e_ol(m,l,k)=dsb(m,k)-waux'*xaux;   % error at subband m
	  sig_ol(m)=(1-a)*sig_ol(m)+a*(xsb(m,k))^2;
	  unlms_ol=2*u/(gamma+(Nw+1)*sig_ol(m));
	  w_ol(m,:)=w_ol(m,:)+unlms_ol*e_ol(m,l,k)*shiftdim(xaux,-1); % new coefficient vector for subband m
      end;
    end;
    wk_ol(l,:,:)=w_ol;
end;

for m=1:M
    mse_ol(m,:)=mean(e_ol(m,:,:).^2,2);    % MSE at subband m
end;
MSE=mean(mse_ol,1);                 % overall MSE
MSEsub = mse_ol;


% Output:
figure,
plot(10*log10(MSE));
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

