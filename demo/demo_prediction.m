% Run a kernel adaptive filtering algorithm for time series prediction.
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

close all;
clear all;

%% PARAMETERS
% Instructions: 1. Uncomment one datafile and one algorithm 2. Execute.

datafile = 'lorenz.dat'; L = 6; N = 10000; horizon = 1;
kaf = aldkrls(struct('nu',1E-4,'M',Inf,'kerneltype','gauss','kernelpar',32));
% kaf = swkrls(struct('c',1E-6,'M',100,'kerneltype','gauss','kernelpar',32));
% kaf = fbkrls(struct('lambda',1E-6,'M',100,'kerneltype','gauss','kernelpar',32));
% kaf = norma(struct('lambda',1E-4,'tau',500,'mu',0.1,'kerneltype','gauss','kernelpar',32));

% datafile = 'mg30.dat'; L = 11; N = 5000; horizon = 1;
% kaf = swkrls(struct('c',1E-6,'M',200,'kerneltype','gauss','kernelpar',0.6));
% kaf = fbkrls(struct('lambda',1E-6,'M',200,'kerneltype','gauss','kernelpar',0.6));
% kaf = norma(struct('lambda',1E-2,'tau',500,'mu',0.5,'kerneltype','gauss','kernelpar',0.6));

%% PROGRAM
tic

data = load(datafile); data = data(:); N = min(N,length(data)-1);
x = data(1:N); X = zeros(N,L);
for i = 1:L, X(i:N,i) = x(1:N-i+1); end % time embedding
Y = data(1+horizon:N+horizon); % desired output

fprintf(1,'Running prediction algorithm')
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

figure; plot(10*log10(SE)); xlabel('samples'); ylabel('squared error (dB)');

figure; hold all; plot(Y); plot(Y_est);
legend('original','prediction'); title(datafile);
