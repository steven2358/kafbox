% Recursive Least-Squares Algorithm with Forgetting Factor
%
% From S. Haykin, "Adaptive Filtering Theory (3rd Ed.)", Prentice Hall,
% Chapter 13.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef rls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        lambda = .99; % forgetting factor
        c = 1E-6; % regularization
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        P = []; % inverse correlation matrix
        w = []; % filter coefficients
    end
    
    methods
        
        function kaf = rls(parameters) % constructor
            if (nargin > 0)
                kaf.lambda = parameters.lambda;
                kaf.c = parameters.c;
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if numel(kaf.w)>0
                y_est = x*kaf.w;
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if numel(kaf.w)==0, % initialize
                m = length(x);
                kaf.w = zeros(m,1);
                kaf.P = kaf.c\eye(m);
            end
            
            k = kaf.P*x'/(kaf.lambda+x*kaf.P*x');
            z = y - x*kaf.w;
            kaf.w = kaf.w + k*z;
            kaf.P = kaf.lambda\(kaf.P - k*x*kaf.P);
        end
        
    end
end
