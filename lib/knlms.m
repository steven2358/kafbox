% Kernel Normalized Least-Mean-Square algorithm with Coherence Criterion
%
% C. Richard, J.C.M. Bermudez, and P. Honeine, "Online Prediction of Time
% Series Data With Kernels," IEEE Transactions on Signal Processing,
% vol. 57, no. 3, pp. 1058-1067, March 2009,
% http://dx.doi.org/10.1109/TSP.2008.2009895
%
% Comment: memories are initialized empty in this implementation
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef knlms < handle
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        eta = .5; % step size
        mu0 = .95; % coherence criterion threshold
        eps = 1E-2; % regularization
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        dict = []; % dictionary
        modict = []; % modulus of the dictionary elements
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = knlms(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if ismember(fn,fieldnames(kaf)),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function train(kaf,x,y) % train the algorithm
            if size(kaf.dict,2)==0 % initialize
                k = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                kaf.dict = x;
                kaf.modict = sqrt(k);
                kaf.alpha = 0;
            else
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                kx = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                C = k./(sqrt(kx)*kaf.modict); % coherence
                if (max(C) <= kaf.mu0), % coherence criterion
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.modict = [kaf.modict; sqrt(kx)];
                    kaf.alpha = [kaf.alpha; 0]; % order increase
                end
            end
            
            h = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            kaf.alpha = kaf.alpha + ...
                kaf.eta / (kaf.eps + h*h') * (y - h*kaf.alpha) * h';
        end
        
    end
end
