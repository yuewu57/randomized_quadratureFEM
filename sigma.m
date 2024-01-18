function [result] =  sigma(MCi)

result=mean(sqrt(MCi(:,1))+0.000001);
%result=mean(abs(sin(MCi(:,1).*MCi(:,2)))+0.2);