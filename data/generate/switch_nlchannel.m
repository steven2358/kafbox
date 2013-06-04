function [x,y] = switch_nlchannel(m,n,p,fun,snr)
% Generate SWITCH_NLCHANNEL data set: input-output data of a nonlinear
% channel whose linear part is changed abruptly each n samples. The
% nonlinear channel is a Wiener system with a randomly chosen linear part.
% 
% Inputs:
% m: number of switches
% n: samples between every switch
% ll: linear channel length
% fun: nonlinear function
% snr: signal-to-noise ratio of additve output noise
%
% Outputs:
% x: system input
% y: system output
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

rs = 1;
randn('state',rs); rand('state',rs); %#ok<RAND>

%% DEFAULT PARAMETERS

if nargin<5, snr = 30; end
if nargin<4, fun = 'tanh(1.5*x)'; end
if nargin<3, p = 5; end
if nargin<2, n = 1000; end
if nargin<1, m = 100; end

%% PROGRAM

f = inline(fun); % Wiener system nonlinearity

N = (m+1)*n;

x = randn(N,1);
X = zeros(N,p);
for i = 1:p,
    X(i:N,i) = x(1:N-i+1);	% time-embedding, each row is a datum
end

xf = zeros(N,1);
for i=1:(m+1),
    sti = (i-1)*n+1;
    ndi = i*n;
    
    ch = [1; randn(p-1,1)]; % linear channel coefficients
    xf(sti:ndi) = X(sti:ndi,:)*ch;
end

% apply nonlinearity
ynn = f(xf); % no noise yet

% add noise
noisevar = 10^(-snr/10)*var(ynn);
noise = sqrt(noisevar)*randn(N,1);

y = ynn + noise;