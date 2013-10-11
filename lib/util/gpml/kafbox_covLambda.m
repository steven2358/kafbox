function [A, B] = kafbox_covLambda(logtheta, xt, z)
% KAFBOX_COVLAMBDA covariance function for an AR(1) process with fixed
% power, parameterized as k(x^(t),x^(t+n)) = lambda^|n|. Uses the
% covariance functions constructors from GPML toolbox version 2.0.
%
% INPUT:	- logtheta = [ log(lambda/(1-lambda))] % hyperparameter
%			- xt: input data, each row is a data point. First column
%			contains temporal indices, other columns are "spatial" data.
%			- z: test set
% OUTPUT:	- A, B: depend on the number of input arguments
% USAGE: [A, B] = kaf_covLambda(logtheta, x, z)
%
% Author: Miguel Lazaro Gredilla, 2012
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

if nargin == 0, A = '1'; return; end              % report number of parameters

t = xt(:,1);           % calculate covariance function only on temporal indices

lambda = exp(logtheta(1))/(1+exp(logtheta));
sf2 = 1;    
n = size(t,1);

if nargin == 2
    A = sf2*lambda.^abs(repmat(t,1,n)-repmat(t',n,1));
elseif nargout == 2                              % compute test set covariances
    z = z(:,1);
    ntst = size(z,1);
    A = sf2*ones(ntst,1);
    B = sf2*lambda.^abs(repmat(t,1,ntst)-repmat(z',n,1));
else                                                % compute derivative matrix
    if z == 1   % wrt lambda
        absn = abs(repmat(t,1,n)-repmat(t',n,1));
        A = sf2*absn.*lambda.^absn.*(1-lambda);
    elseif z == 2   % wrt sf2
        display('It is fixed, you should not be trying to compute this.')
    end
end
