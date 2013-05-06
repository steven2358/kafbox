% Sliding-Window Kernel Recursive Least Squares algorithm
% Author: Steven Van Vaerenbergh, 2013
% Reference: http://dx.doi.org/10.1109/ICASSP.2006.1661394
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef swkrls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        M = 100; % dictionary size
        c = 1E-4; % regularization parameter
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        dicty = []; % output dictionary
        alpha = []; % expansion coefficients
        Kinv = []; % inverse kernel matrix
        prune = false; % flag
    end
    
    methods
        
        function kaf = swkrls(parameters) % constructor
            if (nargin > 0)
                kaf.M = parameters.M;
                kaf.c = parameters.c;
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
            kaf.dict = [kaf.dict; x]; % grow
            kaf.dicty = [kaf.dicty; y];	% grow
            kaf = kaf.grow_kernel_matrix(x); % grow
            
            kaf.prune = false;
            if (size(kaf.dict,1) > kaf.M)
                kaf.prune = true;
                kaf.dict(1,:) = []; % prune
                kaf.dicty(1) = []; % prune
                kaf = kaf.prune_kernel_matrix(); % prune
            end
            
            kaf.alpha = kaf.Kinv*kaf.dicty;
        end
        
    end
    
    methods (Access = 'private')
        
        function kaf = grow_kernel_matrix(kaf,x)
            % calculate inverse of expanded matrix K = [K_inv b;b' d]
            k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
            b = k(1:end-1);
            d = k(end) + kaf.c; % add regularization
            if numel(b)>1
                g_inv = d - b'*kaf.Kinv*b;
                g = 1/g_inv;
                f = -kaf.Kinv*b*g;
                E = kaf.Kinv - kaf.Kinv*b*f';
                kaf.Kinv = [E f;f' g];
            else
                kaf.Kinv = 1/d;
            end
        end
        
        function kaf = prune_kernel_matrix(kaf)
            % calculate inverse of pruned kernel matrix Kp, K = [a b';b Kp]
            m = size(kaf.Kinv,1);
            G = kaf.Kinv(2:m,2:m);
            f = kaf.Kinv(2:m,1);
            e = kaf.Kinv(1,1);
            kaf.Kinv = G - f*f'/e;
        end
        
    end
end
