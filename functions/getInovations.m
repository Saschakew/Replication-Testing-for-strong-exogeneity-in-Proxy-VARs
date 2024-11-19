function [e] = getInovations(u,B)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
e =   u*inv(B)';
end

