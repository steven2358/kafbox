% Fixed-Budget Kernel Recursive Least Squares algorithm
%
% S. Van Vaerenbergh, I. Santamaria, W. Liu, and J.C. Principe, "Fixed-
% budget kernel recursive least-squares," 2010 IEEE International
% Conference on Acoustics Speech and Signal Processing (ICASSP), pp. 1882-
% 1885, 14-19 March 2010, http://dx.doi.org/10.1109/ICASSP.2010.5495350
%
% Comment: label update is not implemented (mu=0)
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef fbkrls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        M = 100; % dictionary size
        lambda = 1E-2; % regularization parameter
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        dicty = []; % output dictionary
        alpha = []; % expansion coefficients
        Kinv = []; % inverse kernel matrix
    end
    
    methods
        
        function kaf = fbkrls(parameters) % constructor
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
            kaf.dict = [kaf.dict; x]; % grow
            kaf.dicty = [kaf.dicty; y];	% grow
            kaf = kaf.grow_kernel_matrix(x); % grow
            
            kaf.alpha = kaf.Kinv*kaf.dicty;
            if (size(kaf.dict,1) > kaf.M)
                ape = abs(kaf.alpha)./diag(kaf.Kinv); % a posteriori error
                [~,ind] = min(ape);
                
                kaf.dict(ind,:) = []; % prune
                kaf.dicty(ind) = []; % prune
                kaf = kaf.prune_kernel_matrix(ind); % prune
            end
            
            kaf.alpha = kaf.Kinv*kaf.dicty;
        end
        
    end
    
    methods (Access = 'private')
        
        function kaf = grow_kernel_matrix(kaf,x)
            % calculate inverse of expanded matrix K = [K_inv b;b' d]
            k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
            b = k(1:end-1);
            d = k(end) + kaf.lambda; % add regularization
            if numel(b)>0
                g_inv = d - b'*kaf.Kinv*b;
                g = 1/g_inv;
                f = -kaf.Kinv*b*g;
                E = kaf.Kinv - kaf.Kinv*b*f';
                kaf.Kinv = [E f;f' g];
            else
                kaf.Kinv = 1/d;
            end
        end
        
        function kaf = prune_kernel_matrix(kaf,ind)
            % calculate inverse of pruned kernel matrix
            m = size(kaf.Kinv,1);
            noind = 1:m;
            noind(ind) = [];
            G = kaf.Kinv(noind,noind);
            f = kaf.Kinv(noind,ind);
            e = kaf.Kinv(ind,ind);
            kaf.Kinv = G - f*f'/e;
        end
    end
end
