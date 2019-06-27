% Kernel Least-Mean-Square algorithm
%
% W. Liu, P.P. Pokharel, and J.C. Principe, "The Kernel Least-Mean-Square
% Algorithm," IEEE Transactions on Signal Processing, vol. 56, no. 2, pp.
% 543-554, Feb. 2008, http://dx.doi.org/10.1109/TSP.2007.907881
%
% Remark: implementation includes a maximum dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef klms < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private')
        eta = .5; % learning rate
        M = 10000; % maximum dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        alpha = []; % expansion coefficients
    end
    
    methods
        function kaf = klms(parameters) % constructor
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
                y_est = zeros(size(x,1),1); % zeros if not initialized
            end
        end
        
        function train(kaf,x,y) % train the algorithm
            if (size(kaf.dict,1)<kaf.M), % avoid infinite growth
                y_est = kaf.evaluate(x);
                err = y - y_est; % instantaneous error
                kaf.dict = [kaf.dict; x]; % add base to dictionary
                kaf.alpha = [kaf.alpha; kaf.eta*err]; % add new coefficient
            end
        end
        
    end
end
