%
% Modulator part of an OFDM transceiver
%
function stime=OFDMmod(sfreq,N,Nactive,Ncp)
stime=[];
for k=1:Nactive:length(sfreq)
    A1=sfreq(k:k+Nactive-1);    % Take a block of QAM symbols
    A=zeros(N,1);               % and put them at right places in 
    A(2:Nactive/2+1)=A1(1:Nactive/2);   % the vector A, before applying
    A(end-Nactive/2+1:end)=A1(Nactive/2+1:end); % IFFT.
    a=ifft(A);              % modulate/convert to time-domain
    a=[a(N-Ncp+1:end); a];    % add CP
    stime=[stime; a];       % concatenate symbols
end
