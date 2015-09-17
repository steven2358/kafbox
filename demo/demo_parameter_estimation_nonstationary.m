% DEMO_PARAMETER_ESTIMATION_NONSTATIONARY Estimation of the parameters of
% the KRLS-T algorithm in a non-statiory setting.
%
% This demo generates data as y = ft(x) where ft is a time-varying
% function. Then it estimates the optimal parameters of the KRLS-T. The
% estimated parameters are: forgetting factor lambda, regularization c and
% Gaussian kernel width. Kernels other than the Gaussian can be used by
% modifying kafbox_parameter_estimation.m.
%
% Author: Steven Van Vaerenbergh, 2013.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

close all
clear

%% PARAMETERS

N = 500; % number of data
sigma = 1; % spatial kernel width
lambda = .995; % forgetting factor (temporal kernel lambda)
c = 1E-2; % AWGN variance

%% PROGRAM
tic

fprintf('\nGenerating time-varying data with spatio-temporal\n');
fprintf('AR(1) covariance...\n');

X = randn(N,1); % 1D input data
time_ind = 0:N-1;
T = toeplitz(time_ind)/2;
Kt = lambda.^T;
K = kernel(X,X,'gauss',sigma);
Kn = c*eye(N);
K = Kt.*K + Kn;
Z = randn(N,1); % GP generator
M = zeros(N,1); % output data mean
Y = chol(K)'*Z + M; % output data

fprintf('Estimating parameters of KRLS-T for these data...\n\n');
[sigma_est,c_est,lambda_est] = kafbox_parameter_estimation(X,Y);

toc
%% OUTPUT

fprintf('\n');
fprintf('        True    Estimated\n');
fprintf('sigma:  %.4f  %.4f\n',sigma,sigma_est)
fprintf('c:      %.4f  %.4f\n',c,c_est)
fprintf('lambda: %.4f  %.4f\n\n',lambda,lambda_est)
fprintf('\n');

% plot data
figure; plot3(time_ind,X(:,1),Y,'+')
xlabel('time'); ylabel('x'); zlabel('y');
