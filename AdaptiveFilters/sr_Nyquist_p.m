%
% Square-root raised-cosine pulse-shape design
%
% h=srNyquist(N,M,alpha);
%
% It designs a square-root raised-cosine pulse-shape with the following
% parameters:
%   N: filter order (filter length = N+1)
%   M: number of samples per symbol period Tb
%   alpha: roll-off factor (between 0 and 1)
%
function h=sr_Nyquist_p(N,M,alpha,gamma);
f2=(1/2/M)*(1+alpha);
h0=sr_cos_p(N,M,alpha);
if rem(N+1,2)==0
    g0=h0(1:(N+1)/2);
else
    g0=h0(1:N/2+1);
end
%
% Set constraint matrices, Sn.
%
Lg=length(g0);
Nconst=1+floor(2*Lg/M);
S=zeros(N+1,N+1,Nconst);
for n=1:Nconst
    for  k=1:N+1
        for l=1:N+1
            if (k-l)==(n-1)*M
                S(k,l,n)=1;
            end
        end
    end
end
%
% Set the matrix Phi
%
Phi=zeros(N+1,N+1);
for k=1:N+1
    for l=1:N+1
        if k==l
            Phi(k,l)=1-2*f2;
        else
            Phi(k,l)=-2*f2*sinc(2*f2*(k-l));
        end
    end
end
%
% Form the matrices Sn' and Phi'.
%
I=eye(Lg);
J=hankel([zeros(Lg-1,1); 1]);
if rem(N+1,2)==1
    J=J(2:end,:);
end
E=[I; J];
Phi1=E'*Phi*E+1e-12;
S1=zeros(Lg,Lg,Nconst);
for n=1:Nconst
    S1(:,:,n)=E'*S(:,:,n)*E;
end
%
% Iterative lease-squares optimization
%
C=chol(Phi1);
g=g0;
for kk=1:100
    B=[];
    for n=1:Nconst
        B=[B; g'*S1(:,:,n)];
    end
    D=[B; gamma*C];
    p=zeros(Nconst+Lg,1);
    p(1)=1;
    g=(g+(D'*D)\(D'*p))/2;
end
h=E*g;






