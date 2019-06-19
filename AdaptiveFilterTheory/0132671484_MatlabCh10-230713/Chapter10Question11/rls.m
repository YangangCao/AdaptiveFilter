function [ W,xp ] = rls( u,d,lambda )
N = min(size(u, 1),size(d, 1));
Nin = size(u,2);
Nout =size(d,2);
w=zeros(Nout,Nin);
W=[];
delta=0.1;
P=eye(Nin)*delta; %reset filter variable between monte carlo runs
for n=1:N
   W=[W;w];
   xp(n,:)=u(n,:)*w';
   e(n,:)=d(n,:)-xp(n,:);     
   kappa=lambda^(-1)*P*u(n,:)'/(1+lambda^(-1)*u(n,:)*P*u(n,:)'); % update kappa as perRLS

   w=w+kappa'*e(n); % update weights

   P=lambda^(-1)*P-lambda^(-1)*u(n,:)*kappa*P; %update as per R

end

