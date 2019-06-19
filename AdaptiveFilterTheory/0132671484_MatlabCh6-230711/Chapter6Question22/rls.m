function [ W,xp ] = rls( u,d,lambda )
N = min(size(u, 1),size(d, 1));
Nin = size(u,2);
Nout =size(d,2);
w=zeros(Nout,Nin);
W=[];
P=eye(Nout)*delta; %reset filter variable between monte carlo runs
for n=1:N
   W=[W;w];
   xp(n,:)=u(n,:)*w';
   e(n,:)=d(n,:)-xp(n,:);     
   kappa=lambda^(-1)*P*xp/(1+lambda^(-1)*xp'*P*xp); % update kappa as perRLS

   W=W+kappa*error(n); % update weights

   P=lambda^(-1)*P-lambda^(-1)*kappa*xp'*P; %update as per R

end

