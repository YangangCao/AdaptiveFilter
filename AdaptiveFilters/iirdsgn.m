%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program let the user to design IIR filters for 		       %
%magnetic recording applications. The target pulse shapes	    %
%are PR4, EPR4 and EEPR4.					                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%To run this program on the student edition of MATLAB, reduce 
%the iteration number "itn" to lower value, such as itn=2000.
%"itn" is defined on line 49.
%
%
% Last updated on April 28, 1998
%

%
%	Design Parameters
%
L=input('over sampling factor (L)? ');	% L=8 is recommended.
D=input('density (D) ? ');
ston=input('signal-to-noise (dB)? ');	%40 dB is usually a good choice.
delay=input('delay? ');	% Typical values are in the range of 1 to 2.
N_n=input('No. of Zeros? ');
N_n1=N_n+1;
N_d=input('No. of poles? ');

%
%	Channel response.
%
h=lorentz(50*L+1,L,D,delay);

h=conv(h,[1 zeros(1,L-1) -1]);
PR=input('PR4 response? (1:PR4, 2:EPR4, 3:EEPR4) ')
if PR==1
	dPR4=[1 0 -1];
elseif PR==2
	dPR4=[1 1 -1 -1];
else
	dPR4=[1 2 0 -2 -1];
end

%
%	FIR equalizer
%
N=32*L+1;
sigman2=sum(h.^2)/(10^(ston/10));
w=leqlzr(h,L,25*L,N,dPR4,sigman2);
c=conv(h,w);

%
%	FIR to IIR conversion
%
itn=10000;
x=filter(h,1,sign(randn(itn,1)));
d=filter(w,1,x);
[A,B]=lsfit(x,d,N_d,N_n1);

%
%	Get the pulse dibit at the equalizer output.
%
f=filter(B,A,h);
figure(1),hold off
c=c(1:length(f));
n_factor=(c*c')/(c*f');
f=f*n_factor;
B=B*n_factor;
n=[1:length(f)];
plot(n,f,'-',n,c,'--')
hold on
plot([1:L:length(c)],round(c(1:L:length(c))),'o')
plot([1 length(f)],[0 0],':')

%
%	residual ISI.
%

efc=f(1:L:length(c))-round(c(1:L:length(c)));
r_isi=efc*efc'/(dPR4*dPR4');
disp(['residual ISI = ' num2str(10*log10(r_isi)) ' dB'])


%
%	Effective S/N at the dector input (equalizer output)
%
S_to_N=10*log10(1/(r_isi+sigman2*(w'*w)/(dPR4*dPR4')));
disp(['Effective S/N at the dector input (equalizer output) = ' num2str(S_to_N) ' dB'])


