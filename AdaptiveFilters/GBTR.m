%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gradient-based timing recovery %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0.9;           % Should be changed for the modified algorithm             
mu=0.0025;          % Step-size parameter
kk=1;xp=0;xm=0;
start=5*L+1;        % To drop the transient of x(t) at the beginning
tau=0.3*ones(1,floor((length(x)-start)/L)); % Initialize the timing offset
            % The timing offset tau is adjusted as the algorithm proceeds.
for k=start:L:length(tau)*L
    tauT=round(tau(kk)*L);
    xp=sqrt(1-beta^2)*x(k+tauT)-beta*xp;
    xm=sqrt(1-beta^2)*x(k+tauT+L/2)-beta*xm;
    tau(kk+1)=tau(kk)+mu*real(xp*xm');
    kk=kk+1;
end