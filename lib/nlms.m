% Normalized Least Mean Squares algorithm
%
% From A. H. Sayed, "Fundamentals of adaptive filtering}", Wiley-IEEE
% Press, 2003, Chapter 5.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef nlms < linear_filter
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu = .9; % step size
        eps = 1E-6; % regularization
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        w = [];
    end
    
    methods
        function obj = nlms(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if ismember(fn,fieldnames(obj)),
                        obj.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(obj,x) % evaluate the algorithm
            if numel(obj.w)>0
                y_est = x*obj.w;
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function train(obj,x,y) % train the algorithm
            if numel(obj.w)==0, % initialize
                obj.w = zeros(length(x),1);
            end
            
            % Algorithm 5.6.1 in reference
            err = (y-x*obj.w);
            obj.w = obj.w + obj.mu/(obj.eps + x*x')*x'*err;
        end
    end
end
