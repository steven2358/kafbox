% Multikernel Normalized Least Mean Square algorithm With Coherence-Based
% Sparsification (MKNLMS-CS)
%
% M. Yukawa, "Multikernel Adaptive Filtering," IEEE Transactions on Signal
% Processing, vol.60, no.9, pp.4672,4682, Sept. 2012,
% http://dx.doi.org/10.1109/TSP.2012.2200889
%
% Comment: version01, August 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef mknlms_cs
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        delta = .9; % coherence criterion threshold
        eta = .1; % step size
        rho = 1E-4; % regularization
        kerneltype = 'gauss'; % kernel type
        kernelpars = [.01 .1 1 10 100 1000]; % kernel parameters
    end
        properties (GetAccess = 'private', SetAccess = 'private') % variables
        dict = []; % dictionary
        alpha = []; % expansion coefficients, H in the original article
        grow = false; % flag
    end
    
    methods
        
        function kaf = mknlms_cs(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if strmatch(fn,fieldnames(kaf),'exact'),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(kaf,X) % evaluate the algorithm
            N = size(X,1);
            if size(kaf.dict,1)>0
                K = multikernel_dict(kaf,X);
                y_est = reshape(K,N,[])*kaf.alpha(:);
            else
                y_est = zeros(N,1);
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            M = length(kaf.kernelpars); % number of distinct kernels
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.alpha = zeros(1,M); % row with coefficients for all kernels
                kaf.grow = true;
            else
                kaf.grow = false;
                K = multikernel_dict(kaf,x);
                if (max(K(:)) <= kaf.delta), % coherence criterion
                    kaf.grow = true;
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.alpha = [kaf.alpha; zeros(1,M)]; % order increase
                end
            end
            K = multikernel_dict(kaf,x);
            kaf.alpha = kaf.alpha + kaf.eta /...
                (kaf.rho + K(:)'*K(:)) * (y - K(:)'*kaf.alpha(:)) * ...
                reshape(K,size(kaf.dict,1),M);
        end
        
        function K = multikernel_dict(kaf,X) % multikernel for dictionary
            M = length(kaf.kernelpars); % number of distinct kernels
            N = size(X,1);
            D = size(kaf.dict,1);
            K = zeros(N,D,M);
            switch kaf.kerneltype
                case 'gauss' % RBF kernel
                    norms1 = sum(X.^2,2);
                    norms2 = sum(kaf.dict.^2,2);
                    mat1 = repmat(norms1,1,D);
                    mat2 = repmat(norms2',N,1);
                    
                    d2 = mat1 + mat2 - 2*X*kaf.dict'; % distance matrix
                    Kalpha = -1./(2*kaf.kernelpars.^2);
                    for m=1:M,
                        K(:,:,m) = exp(d2*Kalpha(m));
                    end
                otherwise	% default case
                    for m=1:M,
                        k = kernel(X,kaf.dict,kaf.kerneltype,kaf.kernelpars(m));
                        K(:,:,m) = k;
                    end
            end
        end
    end
end
