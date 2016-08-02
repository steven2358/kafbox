function [x_embed,y,y_ref,x_test_embed,y_test_ref,H] = ...
    generate_channel_switch(opt)
% Generate CHANNEL_SWITCH data set: input-output data of a nonlinear
% channel whose linear part is changed abruptly at a chosen point. The
% nonlinear channel is a Wiener system with a randomly chosen linear part.
%
% Input options: "opt" is a structure with the following fields:
%   - N: total number of training data points
%   - N_switch: iteration after which the channel switch occurs
%   - N_test: number of test data points (before and after switch)
%   - sigpower: input signal power
%   - chlen: linear channel length
%   - fun: nonlinear function
%   - SNR: signal-to-noise ratio of additve output noise
%
% Outputs:
%   - x_embed: system input with time embedding (each datum is a row)
%   - y: system output
%   - y_ref: noiseless system output
%   - x_test_embed: test system input with time embedding
%   - y_test_ref: noiseless test system output, one column per channel
%   - H: channel matrix (each row is a channel impulse response)
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

%% DEFAULT PARAMETER VALUES

options = struct('N',1500,'N_test',500,'N_switch',500,'chlen',5,...
    'sigpower',1,'fun','tanh(x)','SNR',20);

%% CUSTOM PARAMETER VALUES
if nargin >= 1,
    for opt_name = fieldnames(opt)',
        if strmatch(opt_name,fieldnames(options),'exact'),
            options.(opt_name{1}) = opt.(opt_name{1});
        end
    end
end
    
N = options.N;
N_test = options.N_test;
N_switch = options.N_switch;
sigpower = options.sigpower;
chlen = options.chlen;
fun = options.fun;
SNR = options.SNR;

%% PROGRAM

H = [ones(2,1) randn(2,chlen-1)]; % linear channel coefficients

f = inline(fun); % Wiener system nonlinearity

N_all = N+N_test+chlen-1;
x_all = sqrt(sigpower)*randn(N_all,1);
x_all_embed = zeros(N_all,chlen);
for i = 1:chlen,
    x_all_embed(i:N_all,i) = x_all(1:N_all-i+1); % time-embedding
end
x_all_embed = x_all_embed(chlen:N_all,:);

x_embed = x_all_embed(1:N,:);
x_test_embed = x_all_embed(N+1:N+N_test,:);

% linear filtering
xp1 = x_embed(1:N_switch,:)*H(1,:)';
xp2 = x_embed(N_switch+1:N,:)*H(2,:)';

% apply nonlinearity
y_ref = f([xp1;xp2]); % no noise yet

% add noise
noisevar = 10^(-SNR/10)*var(y_ref);
noise = sqrt(noisevar)*randn(N,1);

y = y_ref + noise;

% get test outputs: 2 columns, no noise
xp3 = x_test_embed*H(1,:)';
xp4 = x_test_embed*H(2,:)';
y_test_ref = f([xp3 xp4]);

% figure;plot([xp1;xp2],y,'.')
