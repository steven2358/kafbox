% Kernel Least-Mean-Square algorithm
%
% W. Liu, P.P. Pokharel, and J.C. Principe, "The Kernel Least-Mean-Square
% Algorithm," IEEE Transactions on Signal Processing, vol. 56, no. 2, pp.
% 543-554, Feb. 2008, http://dx.doi.org/10.1109/TSP.2007.907881
%
% Comment: implementation includes a maximum dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef klms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        eta = .5; % learning rate
        M = 1000; % maximum dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mem = []; % memory
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = klms(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if strmatch(fn,fieldnames(kaf),'exact'),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.mem,1)>0
                k = kernel(kaf.mem,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if (size(kaf.mem,1)<kaf.M), % avoid infinite growth
                y_est = kaf.evaluate(x);
                err = y - y_est;
                kaf.alpha = [kaf.alpha; kaf.eta*err]; % grow
                kaf.mem = [kaf.mem; x]; % grow
            end
        end
        
    end
end
