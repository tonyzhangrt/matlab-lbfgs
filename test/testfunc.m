function [ fval, grad ] = testfunc( x )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
fval = sum(10000*x(1).^4+x(2:end).^4);

grad = 4*x.^3;
grad(1) = 10000*grad(1);

end

