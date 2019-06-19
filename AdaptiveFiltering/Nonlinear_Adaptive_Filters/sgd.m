function [y]=sgd(x,c1,c2)
%  Funcao derivada da sigmoide
%  y = (c2/2*c1)*(c1^2-sgm(x)^2)

y = (c2/2*c1)*(c1^2-sgm(x,c1,c2).^2);
