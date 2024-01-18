function [y] = Fun1(x1,x2)

q=0.49;
eps=10^(-6);


y=abs(x1-x2+eps)^(-q)+10*x2*sign(x2*2-x1);

