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
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        dicty = []; % output dictionary
        alpha = []; % expansion coefficients
        Kinv = []; % inverse kernel matrix
        prune = false; %
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
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            kaf.dict = [kaf.dict; x]; % expand dictionary with row vector
            kaf.dicty = [kaf.dicty; y];	% store subset labels
            
            k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
            m = size(kaf.dict,1);
            kn = k(1:m-1);
            knn = k(end) + kaf.c;                                 % fl{s:1}
            kaf = kaf.inverse_addrowcol(kn,knn); % extend kernel matrix
            
            if (m > kaf.M) % prune dictionary
                kaf.prune = true;
                kaf.dict(1,:) = [];
                kaf.dicty(1) = [];
                kaf = kaf.inverse_removerowcol();
            else
                kaf.prune = false;
            end
            kaf.alpha = kaf.Kinv*kaf.dicty;              %fl{s:n^2-n,m:n^2}
        end
        
        function flops = getflops(kaf) % flops for last iteration
            m = size(kaf.dict,1); % final dictionary size
            n = m - 1;
            floptions = struct(...
                'sum',3*m^2-2*m+2 + kaf.prune*(n^2),...
                'mult',3*m^2+2*m + kaf.prune*(n^2+n),...
                'div',1 + kaf.prune*(1),...
                'kernel',[kaf.kerneltype,m,size(x,2)]);
            flops = kflops(floptions);
        end
        
        function bytes = getbytes(kaf) % bytes used
            m = size(kaf.dict,1);
            bytes = 8*m*(m+2+size(x,2)); % 8 bytes for double precision
        end
    end
    
    methods (Access = 'private')
        
        function kaf = inverse_addrowcol(kaf,b,d)
            % inverse of K = [K_inv b;b' d]
            if numel(b)>1
                g_inv = d - b'*kaf.Kinv*b;     %fl{s:n^2-n+1,n:2*n^2},m=n-1
                g = 1/g_inv;                                       %fl{d:1}
                f = -kaf.Kinv*b*g;                     %fl{s:n^2-n,n:n^2+1}
                E = kaf.Kinv - kaf.Kinv*b*f';              %fl{s:n^2,n:n^2}
                kaf.Kinv = [E f;f' g];
            else
                kaf.Kinv = 1/d;
            end
        end
        
        function kaf = inverse_removerowcol(kaf)
            % inverse of D with K = [a b';b D]
            m = size(kaf.Kinv,1);
            G = kaf.Kinv(2:m,2:m);
            f = kaf.Kinv(2:m,1);
            e = kaf.Kinv(1,1);
            kaf.Kinv = G - f*f'/e;             %fl{s:n^2,n:n^2+n,d:1},n=m-1
        end
    end
end
