function [X,y] = generate_switching_nlchannel(var_options)
% Generate GENERATE_SWITCHING_NLCHANNEL data set: input-output data of a 
% nonlinear channel whose linear part is changed abruptly each n samples. 
% The nonlinear channel is a Wiener system with a randomly chosen linear 
% part.
% 
% Fields of options:
%   - period: samples between every switch
%   - numperiod: number of periods
%   - N: total number of data points (same or less than period*numperiod)
%   - chlen: linear channel length
%   - fun: nonlinear function
%   - SNR: signal-to-noise ratio of additve output noise
%
% Outputs:
% x: system input
% y: system output
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

rs = 1;
randn('state',rs); rand('state',rs); %#ok<RAND>

%% DEFAULT PARAMETER VALUES

options = struct('N',101000,'period',1000,'numperiod',101,'chlen',5,...
    'fun','tanh(1.5*x)','SNR',30);

option_names = fieldnames(options);

%% CUSTOM PARAMETER VALUES
for var_option_name = fieldnames(var_options),
    if any(strmatch(var_option_name,option_names)) % only valid names
        options.var_option_name{1} = var_options.(var_option_name{1});
    end
end

period = options.period;
numperiod = options.numperiod;
N = options.N;
chlen = options.chlen;
fun = options.fun;
SNR = options.SNR;


%% PROGRAM

f = inline(fun); % Wiener system nonlinearity

Np = (numperiod+1)*period;

x = randn(Np,1);
X = zeros(Np,chlen);
for i = 1:chlen,
    X(i:Np,i) = x(1:Np-i+1);	% time-embedding, each row is a datum
end

xf = zeros(Np,1);
for i=1:(numperiod+1),
    sti = (i-1)*period+1;
    ndi = i*period;
    
    ch = [1; randn(chlen-1,1)]; % linear channel coefficients
    xf(sti:ndi) = X(sti:ndi,:)*ch;
end

% apply nonlinearity
ynn = f(xf); % no noise yet

% add noise
noisevar = 10^(-SNR/10)*var(ynn);
noise = sqrt(noisevar)*randn(Np,1);

y = ynn + noise;

N = min(N,Np);
X = X(1:N,:);
y = y(1:N);
