% Kernel Affine Projection algorithm with Coherence Criterion
% Author: Steven Van Vaerenbergh, 2013
% Reference: http://dx.doi.org/10.1109/TSP.2008.2009895
% Comment: memories are initialized empty in this implementation
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef kapcc
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu0 = 1; % coherence criterion threshold
        eta = .1; % step size
        eps = 1E-4; % regularization
        p = 10; % memory length
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mem = []; % memory
        dict = []; % dictionary
        alpha = []; % expansion coefficients
        d = []; % output window
        growmem = false; % flag
        growdict = false; % flag
    end
    
    methods
        
        function kaf = kapcc(parameters) % constructor
            if (nargin > 0)
                kaf.mu0 = parameters.mu0;
                kaf.eta = parameters.eta;
                kaf.eps = parameters.eps;
                kaf.p = parameters.p;
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
            if (length(kaf.d) < kaf.p)
                kaf.mem = [kaf.mem; x]; % grow the memory
                kaf.d = [kaf.d; y]; % grow the memory
                kaf.growmem = true;
            else
                kaf.mem = [kaf.mem(2:end,:); x]; % sliding window
                kaf.d = [kaf.d(2:end); y]; % sliding window
                kaf.growmem = false;
            end
            
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.alpha = 0;
                kaf.growdict = true;
            else
                k = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
                kaf.growdict = false;
                if (max(k) <= kaf.mu0), % coherence criterion
                    kaf.growdict = true;
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.alpha = [kaf.alpha; 0]; % order increase
                end
            end
            
            H = kernel(kaf.mem,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            kaf.alpha = kaf.alpha + ...
                kaf.eta*H'/...
                (kaf.eps*eye(size(H,1)) + H*H')*...
                (kaf.d - H*kaf.alpha);
        end
        
    end
end
