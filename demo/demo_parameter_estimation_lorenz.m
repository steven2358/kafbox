% DEMO_PARAMETER_ESTIMATION_LORENZ Estimation of the parameters of the
% KRLS-T algorithm for predicting the Lorenz attractor time-series.
%
% This demo estimates the optimal parameters of the KRLS-T algorithm when
% applied to predict the Lorenz attractor time-series. The estimated
% parameters are: forgetting factor lambda, regularization c and Gaussian
% kernel width. Kernels other than the Gaussian can be used by modifying
% kafbox_parameter_estimation.m.
%
% Author: Steven Van Vaerenbergh, 2013.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

close all
clear

%% PARAMETERS

horizon = 1; % prediction horizon
embedding = 6; % time-embedding
N = 500; % number of data

%% PROGRAM
tic

fprintf('Loading Lorenz attractor time-series...\n')
[X,Y] = kafbox_data(struct('name','lorenz','horizon',horizon,...
    'embedding',embedding,'N',N));

fprintf('Estimating KRLS-T parameters for %d-step prediction...\n\n',...
    horizon)
[sigma_est,reg_est,lambda_est] = kafbox_parameter_estimation(X,Y);

fprintf('Running KRLS-T with estimated parameters...\n')
Y_est = zeros(N,1);
kaf = krlst(struct('lambda',lambda_est,'M',100,'sn2',reg_est,...
    'kerneltype','gauss','kernelpar',sigma_est));

for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end % progress indicator, 10 dots
    
    Y_est(i) = kaf.evaluate(X(i,:)); % predict the next output
    kaf.train(X(i,:),Y(i)); % train with one input-output pair
    
end
fprintf('\n');
SE = (Y-Y_est).^2; % square error

toc
%% OUTPUT

fprintf('\n');
fprintf('        Estimated\n');
fprintf('sigma:  %.4f\n',sigma_est)
fprintf('c:      %e\n',reg_est)
fprintf('lambda: %.4f\n\n',lambda_est)

fprintf('Average MSE after first 100 samples: %.2fdB\n\n',...
    10*log10(mean(SE(101:end))));

figure; hold all; plot(Y); plot(Y_est);
legend({'original','prediction'},'Location','SE');
title(sprintf('%d-step ahead prediction %s on Lorenz time series',...
    horizon,upper(class(kaf))));
