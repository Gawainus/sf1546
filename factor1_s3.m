function [f1] = factor1(r,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

P = (-r.^3 + 3*r.^2 - 2*r)./(-(1/3)*r.^3 + (3/2)*r.^2 - 2*r + 1)

f1=r./h^2-(1+P)/(2*h);
end
