%   dlcllms.m      
%       Implements the Delayless Closed-Loop Subband (LMS) Adaptive-Filtering Algorithm for REAL valued data.
%       (Algorithm 12.3 - book: Adaptive Filtering: Algorithms and Practical
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
%       Nfd    : Nyquist filter length.
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
Nfd    = 2;


% Design of the fractional delays
% --------------------
K=(Nfd-1)/2;
n=-K:K;
b=sinc(n/L);
win=hamming(Nfd);
h=b.*win';                              % Ideal filter truncated response

% Polyphase decomposition of the Nyquist filter
h=[zeros(1,ceil(Nfd/M)*M-Nfd) h];
Ed=reshape(h,M,length(h)/M);		% analysis filter bank
for m=1:L-1
    Rd(m+1,:)=Ed(m,:);
end
Rd(1,:)=Ed(L,:);
Rd=flipud(Rd);				% synthesis filter bank

F=dftmtx(M);				% DFT matrix
E=(F);
R=(F');


for l=1:Nr
    disp(['Run: ' num2str(l)])
    nc=sqrt(Sn)*randn(1,N*M);		% noise at unknown system output
    xin=sqrt(Sx)*sqrt(12)*(rand(1,N*M)-.5); % input signal
    d=filter(B,A,xin);			% desired signal
    d=d+nc;				% unknown system output

    for m=1:M
	w_cl(m,:)=zeros(1,Nw+1);        % initial coefficient vector for subband m
	x_cl(m,:)=[zeros(1,Nw+1)];
	sig_cl(m)=0;
    end

    d_cl=d;				% desired signal

    Zg=zeros(M*Nw-1,1);
    x_p=zeros(L,1);
    d_p=zeros(L,1);
    xx_frac=zeros(size(Ed,2),M);
    ee_frac=zeros(size(Ed,2),M);
    dd_frac=zeros(size(Ed,2),M);
    yy_frac=zeros(size(Ed,2),M);
    xsb=zeros(M,N);
    for k=1:N
      if mod(k,200)==0 disp(['Iteration: ' num2str(k)]), end;
     % Analysis of input signal
      aux=(k-1)*L+1:k*L;
      x_p=fliplr(xin(aux));
     for m=1:M
	xx_frac(:,m)=cat(1, x_p(m) , xx_frac(1:end-1,m)); % new input signal vector for fractional delay filter of subband m
	x_frac(m)=Ed(m,:)*xx_frac(:,m);	% output of fractional delay filter of subband m
     end;
     xsb(:,k)=E*x_frac';		% input signal split in subbands

     % Mapping of the adaptive filters coefficients into full-band filter
     ww=real(F'*w_cl)/(M);
     G(1,:)=ww(1,1:end-1);
     Dint=(size(Ed,2)-1)/2;
     for m=2:M
	 aux=conv(Ed(m-1,:),ww(m,:));
	 G(m,:)=aux(Dint+2:Dint+1+size(ww,2)-1);
     end;

     GG=reshape(G,1,M*Nw);		% equivalent full-band filter
     [y_dl,Zg]=filter(GG,1,fliplr(x_p),Zg');

      e_aux1=(d_cl((k-1)*L+1:k*L)) - y_dl;
      e_aux1=fliplr(e_aux1);		% L error samples
      for m=1:M
	ee_frac(:,m)=cat(1, (e_aux1(m)) , ee_frac(1:end-1,m));
	e_frac(m)=Ed(m,:)*ee_frac(:,m);
      end;
      esb(:,k)=E*e_frac';		% error samples split in subbands
      for m=1:M				% adaptation in the subbands
	  x_cl(m,:)=[xsb(m,k) x_cl(m,1:Nw)]; % new input signal vector for subband m
	  sig_cl(m)=(1-a)*sig_cl(m)+a*abs(xsb(m,k))^2; 
	  unlms_cl=u/(gamma+(Nw+1)*sig_cl(m));
	  w_cl(m,:)=w_cl(m,:)+2*unlms_cl*esb(m,k)'*(x_cl(m,:)); % new coefficient vector for subband m
      end;
      e_cl(l,(k-1)*L+1:k*L)=e_aux1;
    end;
end;


MSE=mean(abs(e_cl).^2,1);		% MSE


% Output:
figure,
plot(10*log10(MSE));
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

