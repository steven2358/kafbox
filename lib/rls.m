% Recursive Least-Squares Algorithm with exponential weighting
%
% From S. Haykin, "Adaptive Filtering Theory (3rd Ed.)", Prentice Hall,
% Chapter 13.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef rls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        lambda = .99; % forgetting factor
        c = 1E-4; % regularization
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        P = []; % inverse correlation matrix
        w = []; % filter coefficients
    end
    
    methods
        
        function obj = rls(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if ismember(fn,fieldnames(kaf)),
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
        
        function obj = train(obj,x,y) % train the algorithm
            if numel(obj.w)==0, % initialize
                m = length(x);
                obj.w = zeros(m,1);
                obj.P = obj.c\eye(m);
            end
            
            k = obj.P*x'/(obj.lambda+x*obj.P*x');
            z = y - x*obj.w;
            obj.w = obj.w + k*z;
            obj.P = obj.lambda\(obj.P - k*x*obj.P);
        end
        
    end
end
