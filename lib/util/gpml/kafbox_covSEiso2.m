function [A, B] = kafbox_covSEiso2(loghyper, xt, z)
% KAFBOX_COVSEISO2 Squared Exponential covariance function with isotropic
% distance measure. The covariance function is parameterized as:
% k(x^p,x^q) = sf2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2), where the P 
% matrix is ell^2 times the unit matrix and sf2 is the signal variance.
%
% INPUT:	- loghyper = [log(ell); log(sqrt(sf2))] % hyperparameters
%			- xt: input data, each row is a data point. First column
%			contains temporal indices, other columns are "spatial" data.
%			- z: test set
% OUTPUT:	- A, B: depend on the number of input arguments
% USAGE: [A, B] = kaf_covSEiso2(loghyper, xt, z)
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/
%
% Original covSEiso.m (C) Copyright 2006 by Carl Edward Rasmussen as part
% of the GPML toolbox version 2.0. The code and associated documentation
% are avaiable from http://www.GaussianProcess.org/gpml/code

if nargin == 0, A = '2'; return; end              % report number of parameters

x = xt(:,2:end);           % calculate covariance function only on spatial data

ell = exp(loghyper(1));                           % characteristic length scale
sf2 = exp(2*loghyper(2));                                     % signal variance

if nargin == 2
  A = sf2*exp(-kafbox_sq_dist(x'/ell)/2);
elseif nargout == 2                              % compute test set covariances
  z = z(:,2:end);
  A = sf2*ones(size(z,1),1);
  B = sf2*exp(-kafbox_sq_dist(x'/ell,z'/ell)/2);
else                                                % compute derivative matrix
  if z == 1                                                   % first parameter
    A = sf2*exp(-kafbox_sq_dist(x'/ell)/2).*kafbox_sq_dist(x'/ell);  
  else                                                       % second parameter
    A = 2*sf2*exp(-kafbox_sq_dist(x'/ell)/2);
  end
end
