function [y]=sgm(x,c1,c2)
%  Funcao sigmoide
%  y = 2*c1/(1+exp(-c2*x))-c1

y=2*c1./(1+exp(-c2*x))-c1;
