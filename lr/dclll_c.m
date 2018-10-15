function [ H_t T ] = dclll_c( H, delta )
%CLLC_C Summary of this function goes here
%   A wrapper method to call the c implementation of CLLL

if (nargin <= 1)
    delta = 0.75;
end

H_D = pinv(H');

[Q R] = qr(H_D, 0);

T = eye(size(R, 2));
[T] = clll_core_c(R, T, delta);
T = inv(T)';
H_t = H * T;
end