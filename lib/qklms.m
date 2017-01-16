% Quantized Kernel Least Mean Square algorithm
%
% B. Chen, S. Zhao, P. Zhu, and J.C. Principe, "Quantized Kernel Least Mean
% Square Algorithm," IEEE Transactions on Neural Networks and Learning
% Systems, vol. 23, no. 1, pp. 22-32, Jan. 2012,
% http://dx.doi.org/10.1109/TNNLS.2011.2178446
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef qklms < handle
    
    properties (GetAccess = 'public', SetAccess = 'private')
        eta = .9; % learning rate
        epsu = .1;
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % codebook (dictionary)
        alpha = []; % expansion coefficients
    end
    
    methods
        function kaf = qklms(parameters) % constructor
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
            y_est = kaf.evaluate(x); % evaluate function output
            err = y - y_est; % instantaneous error
            
            m = size(kaf.dict,1);
            if m==0
                d2 = kaf.epsu^2 + 1; % force addition of initial base
            else % find distance to closest dictionary element
                [d2,j] = min(sum((kaf.dict - repmat(x,m,1)).^2,2));
            end
            if d2 <= kaf.epsu^2, % reduced coefficient update
                kaf.alpha(j) = kaf.alpha(j) + kaf.eta*err;
            else
                kaf.dict = [kaf.dict; x]; % add base to dictionary
                kaf.alpha = [kaf.alpha; kaf.eta*err]; % add new coefficient
            end
        end
    end
end
