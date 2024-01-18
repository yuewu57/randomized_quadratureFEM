function [y] = Fun_rough(x1,x2)

q=0.49;


y=abs(x1-x2)^(-q)+10*sin(2^3*pi*x1)*sign(x2*2-x1);
