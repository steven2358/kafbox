% Kernel Least Mean Squares algorithm with Coherence-Sparsification
% criterion and Adaptive L1-norm regularization
%
% Wei Gao, Jie Chen, Cedric Richard, Jianguo Huang, and Remi Flamary,
% "Kernel LMS algorithm with forward-backward splitting for dictionary
% learning,", 2013 IEEE International Conference on Acoustics, Speech, and
% Signal Processing (ICASSP 2013), Vancouver, Canada, March 2013.
% http://dx.doi.org/10.1109/ICASSP.2013.6638763
%
% Comment: Code contributed by Wei Gao and Cedric Richard.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef klms_csal1
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        eta = .1; % step-size
        mu0 = .95; % threshold for coherence criterion
        lambda = 5E-4; % sparsification threshold
        eps_alpha = 1E-6; % constant to prevent denominator from vanishing
        kerneltype = 'gauss'; % kernel type
        kernelpar = .5; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        dict = []; % dictionary
        modict = []; % modulus of the dictionary elements
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = klms_csal1(parameters) % constructor
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
                kd = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                kaf.dict = x;
                kaf.modict = sqrt(kd);
                kaf.alpha = y/kd;
                w = 1./abs(kaf.alpha + kaf.eps_alpha); % for sparsification
            else
                w = 1./abs(kaf.alpha + kaf.eps_alpha); % for sparsification
                
                kd = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                kx = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                C = kd./(sqrt(kx)*kaf.modict); % coherence
                
                if (max(C) <= kaf.mu0), % coherence criterion
                    kaf.dict = [kaf.dict; x]; % order increase
                    kaf.modict = [kaf.modict; sqrt(kx)];
                    kaf.alpha = [kaf.alpha; 0]; % order increase
                    kx = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                    kd = [kd; kx];
                    w = [w; 1];
                end
            end
            
            y_est = kaf.evaluate(x);
            e = y - y_est;
            kaf.alpha = kaf.alpha + kaf.eta * kd * e;
            
            kaf.alpha = kaf.prox(kaf.lambda,kaf.alpha,w);
            
            % prune
            idx = find(kaf.alpha == 0);
            if ~isempty(idx)
                kaf.dict(idx,:) = [];
                kaf.modict(idx,:) = [];
                kaf.alpha(idx,:) = [];
            end
        end
        
    end
    
    methods (Static = true)
        
        % proximal operator for l1 norm
        function alphap = prox(lambda,alpha,w)
            alphap = sign(alpha).*max(abs(alpha) - lambda*w, 0);
        end
        
    end
end
