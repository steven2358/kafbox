% Naive Online regularized Risk Minimization Algorithm
%
% J. Kivinen, A.J. Smola, and R.C. Williamson, "Online learning with
% kernels," IEEE Transactions on Signal Processing, vol. 52, no. 8,
% pp. 2165-2176, Aug. 2004, http://dx.doi.org/10.1109/TSP.2004.830991
%
% Remark: using squared loss function
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef norma < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private')
        tau = 500; % memory size (terms retained in truncation)
        lambda = 1E-2; % regularization parameter
        eta = .5; % learning rate
        tcoeff = 0; % learning rate coefficient: eta_t = eta * t^tcoeff
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        alpha = []; % expansion coefficients
    end
    
    methods        
        function kaf = norma(parameters) % constructor
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
            kaf.alpha = (1-kaf.lambda*kaf.eta)*kaf.alpha;
            
            y_est = kaf.evaluate(x);
            err = y - y_est;
            
            kaf.alpha = [kaf.alpha; kaf.eta*err]; % grow
            kaf.dict = [kaf.dict; x]; % grow
            if length(kaf.alpha)>kaf.tau
                kaf.alpha(1) = []; % prune
                kaf.dict(1,:) = []; % prune
            end
        end
    end
end
