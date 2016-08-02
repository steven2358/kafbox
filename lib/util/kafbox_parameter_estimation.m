function [sigma,reg,lambda] = ...
    kafbox_parameter_estimation(X,Y,time_ind,epochs)
% KAFBOX_PARAMETER_ESTIMATION Gaussian-process based estimation of the
% parameters of the KRLS-T algorithm.
%
% Steven Van Vaerenbergh, Ignacio Santamaria, and Miguel Lazaro-Gredilla,
% "Estimation of the forgetting factor in kernel recursive least squares."
% 2012 IEEE International Workshop on Machine Learning for Signal
% Processing (MLSP), 2012. http://dx.doi.org/10.1109/MLSP.2012.6349749
%
% INPUT:	- X: input data, each row is a data point
%			- Y: output data, one column.
%           - time_ind: temporal indices of points
%           - epochs: number of epochs for minimization
% OUTPUT:	- sigma: estimated kernel length scale (width)
%           - reg: estimated regularization
%           - lambda: estimated forgetting factor.
% USAGE: [sigma,reg,lambda] =
%             kafbox_parameter_estimation(X,Y,time_ind,epochs)
%
% Author: Steven Van Vaerenbergh, 2013.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

N = size(X,1);
if nargin<4
    epochs = 15;
end
if nargin<3
    time_ind = (0:N-1)'; % assume that samples are taken on synchronous time instants
end
Xt = [time_ind X];

% initialization values
noisevar0 = mean(Y)^2/4;
lambda0 = 0.99;
ell0 = 1;
sf20 = 1;

% The standard covariance function in kernel adaptive filtering is a
% Gaussian kernel with additive noise. Forgetting is introduced through an 
% AR(1) process with forgetting factor lambda. See publication for details.
covfunc = {'kafbox_covSum',{{'kafbox_covProd',...
    {'kafbox_covLambda', 'kafbox_covSEiso2'}},'kafbox_covNoise'}};
loghyper0 = [log(lambda0/(1-lambda0)); log(ell0); log(sqrt(sf20)); 0.5*log(noisevar0)];

% Minimize negative log maximal likelihood using GPML toolbox
loghyper1 = kafbox_minimize(loghyper0,'kafbox_gpr',epochs,covfunc,Xt,Y);

lambda_est = exp(loghyper1(1))/(1+exp(loghyper1(1)));
lambda_est = lambda_est^2;	% correspondance between st-gp and krls-t
ell_est = exp(loghyper1(2));
sf2_est = exp(2*loghyper1(3));
noisevar_est = exp(2*loghyper1(4));
 
sigma = ell_est;	% length-scale
reg = noisevar_est/sf2_est;  % noise-to-signal parameter (regularization)
lambda = lambda_est;   % forgetting factor
