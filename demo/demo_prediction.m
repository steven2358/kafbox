% Run a kernel adaptive filtering algorithm for time series prediction.
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/
clear all;

%% PARAMETERS

swkrls_pars = struct('c',1E-6,'M',100,'kerneltype','gauss','kernelpar',32);
datafile = 'lorenz.dat'; L = 6; N = 10000; horizon = 1;

%% PROGRAM
tic

data = load(datafile); data = data(:);
x = data(1:N); X = zeros(N,L);
for i = 1:L, X(i:N,i) = x(1:N-i+1); end % time embedding
Y = data(1+horizon:N+horizon); % desired output

fprintf(1,'Running prediction algorithm...\n')
kaf = swkrls(swkrls_pars);
Y_est = zeros(N,1);
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end
    
    Y_est(i) = kaf.evaluate(X(i,:)); % make prediction
    kaf = kaf.train(X(i,:),Y(i)); % train
end
fprintf('\n');

SE = (Y-Y_est).^2; % test error

toc
%% OUTPUT

fprintf('MSE after first 1000: %.2fdB\n\n',10*log10(mean(SE(1001:end))));

figure; hold all; plot(Y); plot(Y_est);