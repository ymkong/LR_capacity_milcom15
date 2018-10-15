function [ H_t, T ] = sa_c( H )
%SA_C Summary of this function goes here
%   A wrapper of SA's algorithm (C implementation)
%   Written By Qi, 5/8/2011

G = H' * H;

C = inv(G);

T = eye(size(C, 1));

T = sa_core_c(G, C, T);

H_t = H * T;

end

