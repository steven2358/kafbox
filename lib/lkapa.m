% Leaky Kernel Affine Projection Algorithm
%
% W. Liu and J.C. Principe, "Kernel Affine Projection Algorithms", EURASIP
% Journal on Advances in Signal Processing, Volume 2008, Article ID 784292,
% 12 pages. http://dx.doi.org/10.1155/2008/784292
%
% Comment: This implementation includes a maximum dictionary size M. With
% M=Inf this algorithm is equivalent to KAPA-3 from the publication. With
% M=Inf and lambda=0 it is equivalent to KAPA-1.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef lkapa
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        eta = .05; % learning rate
        lambda = 1E-2; % Tikhonov regularization
        M = 1000; % maximum dictionary size
        P = 20; % number of regressors
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        xmem = []; % input memory
        ymem = []; % output memory
        dict = []; % dictionary
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = lkapa(parameters) % constructor
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
                y_est = zeros(size(x,1),1);
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.alpha = kaf.eta*y;
                kaf.xmem = x;
                kaf.ymem = y;
            else
                if size(kaf.dict,1) < kaf.M,
                    if size(kaf.xmem,1)<kaf.P
                        % grow memory
                        kaf.xmem = [kaf.xmem; x];
                        kaf.ymem = [kaf.ymem; y];
                    else
                        % slide memory
                        kaf.xmem = [kaf.xmem(2:kaf.P,:); x];
                        kaf.ymem = [kaf.ymem(2:kaf.P); y];
                    end
                    
                    ymem_est = kaf.evaluate(kaf.xmem);
                    e = kaf.ymem - ymem_est;
                    
                    kaf.dict = [kaf.dict; x]; % grow dictionary
                    kaf.alpha = [(1-kaf.lambda*kaf.eta)*kaf.alpha; 0];
                    % leak and grow coefficients
                    
                    m = size(kaf.alpha,1);
                    p = size(kaf.xmem,1);
                    
                    % update p last coefficients
                    kaf.alpha(m-p+1:m) = kaf.alpha(m-p+1:m) + kaf.eta*e;
                end
            end
        end
        
    end
end
