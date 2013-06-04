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
% http://sourceforge.net/projects/kafbox/

classdef knlms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu0 = .9; % coherence criterion threshold
        eta = .1; % step size
        eps = 1E-4; % regularization
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        alpha = []; % expansion coefficients
        grow = false; % flag
    end
    
    methods
        
        function kaf = knlms(parameters) % constructor
            if (nargin > 0)
                kaf.mu0 = parameters.mu0;
                kaf.eta = parameters.eta;
                kaf.eps = parameters.eps;
                kaf.kerneltype = parameters.kerneltype;
                kaf.kernelpar = parameters.kernelpar;
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.alpha = 0;
                kaf.grow = true;
            else
                k = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
                kaf.grow = false;
                if (max(k) <= kaf.mu0), % coherence criterion
                    kaf.grow = true;
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.alpha = [kaf.alpha; 0]; % order increase
                end
            end
            
            h = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            kaf.alpha = kaf.alpha + ...
                kaf.eta / (kaf.eps + h*h') * (y - h*kaf.alpha) * h';
        end
        
    end
end
