% Normalized Least Mean Squares algorithm
%
% From A. H. Sayed, "Fundamentals of adaptive filtering}", Wiley-IEEE Press,
% 2003, Chapter 5.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef nlms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu = .9; % step size
        eps = 1E-6; % regularization
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        w = [];
    end
    
    methods
        
        function obj = nlms(parameters) % constructor
            if (nargin > 0)
                obj.mu = parameters.mu;
                obj.eps = parameters.eps;
            end
        end
        
        function y_est = evaluate(obj,x) % evaluate the algorithm
            if numel(obj.w)>0
                y_est = x*obj.w;
            else
                y_est = 0;
            end
        end
        
        function obj = train(obj,x,y) % train the algorithm
            if numel(obj.w)==0, % initialize
                obj.w = zeros(length(x),1);
            end
            
            % Eq. (5.6.4) in reference
            obj.w = obj.w + obj.mu/(obj.eps + x*x')*x'*(y-x*obj.w);
        end
        
    end
end
