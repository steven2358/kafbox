% Kernel Affine Projection (KAP) algorithm with coherence criterion
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

classdef kap
    
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
        memd = []; % output memory
        dict = []; % dictionary
        alpha = []; % expansion coefficients
        growmem = false; % flag
        growdict = false; % flag
    end
    
    methods
        
        function kaf = kap(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if strmatch(fn,fieldnames(kaf),'exact'),
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
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if (length(kaf.memd) < kaf.p)
                kaf.memx = [kaf.memx; x]; % grow the memory
                kaf.memd = [kaf.memd; y]; % grow the memory
            else
                kaf.memx = [kaf.memx(2:end,:); x]; % sliding memory
                kaf.memd = [kaf.memd(2:end); y]; % sliding memory
            end
            
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.alpha = 0;
            else
                k = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
                if (max(k) <= kaf.mu0), % coherence criterion
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.alpha = [kaf.alpha; 0]; % order increase
                end
            end
            
            H = kernel(kaf.memx,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            kaf.alpha = kaf.alpha + ...
                kaf.eta*H'/...
                (kaf.eps*eye(size(H,1)) + H*H')*...
                (kaf.memd - H*kaf.alpha);
        end
        
    end
end
