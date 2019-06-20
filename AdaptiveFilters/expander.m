%
% EXPANDER: y=expander(x,L) 
% When x is a vector, this function adds L-1 after each element of x. 
% When x is a matrix, each column of it treated as vector and expanded L
% fold.
function y=expander(x,L)
[M,N]=size(x);
if (N==1)|(M==1)
    if M<N
        y=zeros(1,N*L);
        y(1:L:end)=x;
    else
        y=zeros(M*L,1);
        y(1:L:end)=x;
    end
else
    y=zeros(M*L,N);
    y(1:L:end,:)=x;
end
