% Sliding-Window Kernel Recursive Least Squares algorithm
% Author: Steven Van Vaerenbergh, 2013
% Reference: http://dx.doi.org/10.1109/ICASSP.2006.1661394
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef swkrls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        M = 100; % memory size
        c = 1E-4; % regularization parameter
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        dict = []; % dictionary
        dicty = []; % output dictionary
        alpha = []; % expansion coefficients
        Kinv = []; % inverse kernel matrix
    end
    
    methods
        
        function kaf = swkrls(parameters) % constructor
            if(nargin > 0)
                kaf.M = parameters.M;
                kaf.c = parameters.c;
                kaf.kerneltype = parameters.kerneltype;
                kaf.kernelpar = parameters.kernelpar;
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate method
            if size(kaf.dict,1)>0
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train method
            kaf.dict = [kaf.dict; x]; % expand dictionary with row vector
            kaf.dicty = [kaf.dicty; y];	% store subset labels
            
            k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
            m = size(kaf.dict,1);
            kn = k(1:m-1);
            knn = k(end) + kaf.c;
            kaf.Kinv = kaf.inverse_addrowcol(kn,knn,kaf.Kinv); % extend kernel matrix
            
            if (m>kaf.M)
                % prune dictionary
                kaf.dict(1,:) = [];
                kaf.dicty(1) = [];
                kaf.Kinv = kaf.inverse_removerowcol(kaf.Kinv);
            end
            kaf.alpha = kaf.Kinv*kaf.dicty;
        end
    end
    
    methods (Static = true, Access = private)
        
        function K_inv = inverse_addrowcol(b,d,A_inv)
            % returns the inverse of K = [A b;b' d]
            if numel(b)>1
                g_inv = d - b'*A_inv*b;
                g = 1/g_inv;
                f = -A_inv*b*g;
                E = A_inv - A_inv*b*f';
                K_inv = [E f;f' g];
            else
                K_inv = 1/d;
            end
        end
        
        function D_inv = inverse_removerowcol(K_inv)
            % calculates the inverse of D with K = [a b';b D]
            m = size(K_inv,1);
            G = K_inv(2:m,2:m);
            f = K_inv(2:m,1);
            e = K_inv(1,1);
            D_inv = G - f*f'/e;
        end
    end
end