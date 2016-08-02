function [u,d,dref] = generate_richardbench(N)
% Generate RICHARDBENCH signal.
%
% Benchmark signal introduced in C. Richard, J.C.M. Bermudez, P. Honeine,
% "Online Prediction of Time Series Data With Kernels," IEEE Transactions
% on Signal Processing, vol.57, no.3, pp.1058,1067, March 2009.
% 
% Comment: copyright Cedric Richard, http://cedric-richard.fr/
%
% Input: N: number of data points
% 
% Outputs: u: input sequence (1-dimensional sequence [u(:,1)])
%          d: noisy desired output (1-dimensional sequence)
%          dref: noise-free desired output
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

v = zeros(1,N+1);
v(1) = 0.5;
u = 0.25*randn(N+1,1);
dref = zeros(1,N+1);

for t=2:N+1,
    v(t) = 1.1*exp(-abs(v(t-1)))+u(t);
    dref(t) = v(t)^2;
end
d = dref+randn(1,N+1);

d(1) = [];
dref(1) = [];
u(1) = [];
 