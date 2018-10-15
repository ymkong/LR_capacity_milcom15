function [ H_t, T, info ] = clll_c( H, info )
%CLLC_C Summary of this function goes here
%   A wrapper method to call the c implementation of CLLL

if (nargin <= 1)
    info = struct('delta', 0.75);
end

[Q, R] = qr(H, 0);

T = eye(size(R, 2));
[T, info] = clll_core_c(R, T, info);
H_t = H * T;

if (nargout >= 2)
    info.R = R;
end

end

