% Multikernel Normalized Least Mean Square algorithm With Coherence-Based
% Sparsification (MKNLMS-CS)
%
% M. Yukawa, "Multikernel Adaptive Filtering," IEEE Transactions on Signal
% Processing, vol.60, no.9, pp.4672,4682, Sept. 2012,
% http://dx.doi.org/10.1109/TSP.2012.2200889
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef mknlmscs
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu0 = .9; % coherence criterion threshold
        eta = .1; % step size
        rho = 1E-4; % regularization
        kerneltype = 'gauss'; % kernel type
        kernelpars = [.01 .1 1 10 100 1000]; % kernel parameters
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        h = []; % expansion coefficients
        grow = false; % flag
    end
    
    methods
        
        function kaf = mknlmscs(parameters) % constructor
            if (nargin > 0)
                kaf.mu0 = parameters.mu0;
                kaf.eta = parameters.eta;
                kaf.rho = parameters.rho;
                kaf.kerneltype = parameters.kerneltype;
                kaf.kernelpars = parameters.kernelpars;
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                K = multikernel_dict(kaf,x);
                y_est = K(:)'*kaf.h(:);
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            M = length(kaf.kernelpars); % number of distinct kernels
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.h = zeros(1,M); % coefficients of all kernels as row (vs. column as in publication)
                kaf.grow = true;
            else
                kaf.grow = false;
                K = multikernel_dict(kaf,x);
                if (max(K(:)) <= kaf.mu0), % coherence criterion
                    kaf.grow = true;
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.h = [kaf.h; zeros(1,M)]; % order increase
                end
            end
            
            K = multikernel_dict(kaf,x);
            kaf.h = kaf.h + ...
                kaf.eta / (kaf.rho + K(:)'*K(:)) * (y - K(:)'*kaf.h(:)) * K;
        end
        
        function K = multikernel_dict(kaf,x) % multikernel for dictionary
            M = length(kaf.kernelpars); % number of distinct kernels
            K = zeros(size(kaf.dict,1),M);
            for m=1:M,
                K(:,m) = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpars(m));
            end
        end
    end
end