%   cfdlms.m      
%       Implements the Constrained Frequency-Domain (LMS) Algorithm for REAL valued data.
%       (Algorithm 12.4 - book: Adaptive Filtering: Algorithms and Practical
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
L      = M/2;                % interpolation/decimation factor
a      = 0.01;
Nw     = 5;
%  M=2*L;                    % number of bins


% Constrained Frequency-Domain LMS
for l=1:Nr
    disp(['Run: ' num2str(l)])
    nc=sqrt(Sn)*randn(1,N*M);     % noise at system output 
    xin=rand(1,N*M);
    xin=sqrt(12)*(xin-.5);        % input signal
    d=filter(B,A,xin);            % desired signal
    d=d+nc;                       % unknown system output

    y=zeros(M,N);
    xsb=zeros(M,N);
    dsb=zeros(L,N);
    Zy=zeros((L)*(Nw+1)-1,1);
    sig=zeros(M,1);
    ww=zeros(M,(Nw+1));
    wwc=zeros((Nw+1),(M));
    w=reshape(ww,Nw+1,M);
    uu=zeros(M,Nw+1);
    xin=[zeros(1,L) xin];
    for k=1:N
      % Serial-to-parallel converter
      aux1=[1:L]+(k-1)*L;
      aux2=[L+1:M]+(k-1)*L;
      x_p=[fliplr(xin(aux2)) fliplr(xin(aux1))]; xsb(:,k)=x_p';
      aux=[1:L]+(k-1)*L;
      d_p=fliplr(d(aux)); dsb(:,k)=d_p';
      ui=fft(xsb(:,k))/sqrt(M);
      uu=[ui , uu(:,1:end-1)];
      for m=1:M
          uy(m,:)=uu(m,:)*ww(m,:).';
      end;
      y(:,k)=ifft(uy)*sqrt(M);
	     
      e=(dsb(:,k))-(y(1:L,k));		% error samples
      out(aux)=y(1:L,k);		% adaptive filter output
      IL=[eye(L) ; zeros(M-L,L)];
      et=fft([e;zeros(M-L,1)])/sqrt(M); % F*IL*e;
      clear wwc;
      for m=1:M
	  sig(m)=(1-a)*sig(m)+a*abs(ui(m)).^2;
	  wwc(m,:)= u/(gamma+(Nw+1)*sig(m))*conj(uu(m,:))*(et(m)); % new increment for the coefficient vector for subband m (unconstrained)
      end;
      % Constraint
      waux=fft(wwc)/sqrt(M);
      wwc=sqrt(M)*ifft([waux(1:L,:); zeros(M-L,Nw+1)]); % new increment for the coefficient vector for subband m (constrained)
      ww=ww+wwc; 	% new coefficient vector for subband m (constrained)
      e_fda(:,k)=abs(e).^2;
    end;
    e_fd(l,:)=e_fda(:)';
end;


MSE=mean(abs(e_fd),1);
MSE=shiftdim(MSE,1);


% Output:
figure,
plot(10*log10(MSE));
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

