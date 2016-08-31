function [ fval, grad ] = rosenbrockfunc( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% N = size(x,1);
% hN = floor(N/2);
a = 1;
b = 100;

% fval = sum((a-x(1:2:1+2*hN)).^2)+ sum (b*()) 
fval = (a-x(1))^2+b*(x(2)-x(1)^2)^2;

grad = [2*(x(1)-a) - 4*b*(x(2)-x(1))*x(1) ; 2*b*(x(2)-x(1)) ]; 



end

