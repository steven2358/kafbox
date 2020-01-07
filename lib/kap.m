% Kernel Affine Projection (KAP) algorithm with coherence criterion
%
% C. Richard, J.C.M. Bermudez, and P. Honeine, "Online Prediction of Time
% Series Data With Kernels," IEEE Transactions on Signal Processing,
% vol. 57, no. 3, pp. 1058-1067, March 2009,
% http://dx.doi.org/10.1109/TSP.2008.2009895
%
% Remark: memories are initialized empty in this implementation
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef kap < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu0 = .95; % coherence criterion threshold
        eta = .5; % step size
        eps = 1E-2; % regularization
        p = 20; % memory length
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        memx = []; % input memory
        memy = []; % output memory
        dict = []; % dictionary
        modict = []; % modulus of the dictionary elements
        alpha = []; % expansion coefficients
    end
    
    methods
        function kaf = kap(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)'
                    if ismember(fn,fieldnames(kaf))
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
            if (length(kaf.memy) < kaf.p)
                kaf.memx = [kaf.memx; x]; % grow the memory
                kaf.memy = [kaf.memy; y]; % grow the memory
            else
                kaf.memx = [kaf.memx(2:end,:); x]; % sliding memory
                kaf.memy = [kaf.memy(2:end); y]; % sliding memory
            end
            
            if size(kaf.dict,2)==0 % initialize
                k = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                kaf.dict = x;
                kaf.modict = sqrt(k);
                kaf.alpha = 0;
            else
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                kx = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                C = k./(sqrt(kx)*kaf.modict); % coherence
                if (max(C) <= kaf.mu0) % coherence criterion
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.alpha = [kaf.alpha; 0]; % order increase
                end
            end
            
            H = kernel(kaf.memx,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            kaf.alpha = kaf.alpha + ...
                kaf.eta*H'/...
                (kaf.eps*eye(size(H,1)) + H*H')*...
                (kaf.memy - H*kaf.alpha);
        end
    end
end
