function [v,d,dref] = generate_doddbench(N)
% Generate DODDBENCH signal.
%
% Benchmark signal introduced in Dodd, T.J., Kadirkamanathan, V. and
% Harrison, R.F., "Function estimation in Hilbert space using sequential
% projections," Proc. of the IFAC Conf. on Intelligent Control Systems and
% Signal Processing, 113-118, 2003.
% 
% Comment: copyright Cedric Richard, http://cedric-richard.fr/
%
% Input: N: number of data points
% 
% Outputs: v: input sequence (2-dimensional sequence [v(:,1);v(:,2)])
%          d: noisy desired output (1-dimensional sequence)
%          dref: noise-free desired output
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

dref = zeros(1,N+2);
dref(1:2)=[0.1 0.1];

for t=3:N+2
    dref(t) = (0.8-0.5*exp(-dref(t-1)^2))*dref(t-1) - ...
        (0.3+0.9*exp(-dref(t-1)^2))*dref(t-2)+0.1*sin(pi*dref(t-1));
end
d = dref + 0.1*randn(1,N+2);
v = [d(1:N); d(2:N+1)]';

d(1:2)=[];
dref(1:2)=[];
