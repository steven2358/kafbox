% Dummy prediction method that always returns the last seen output.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef memory_cell < handle
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        y_mem;
    end
    
    methods
        % dummy constructor
        function obj = memory_cell(~)
        end
        
        % evaluate the algorithm
        function y_est = evaluate(obj,x)
            if ~isempty(obj.y_mem)
                y_est = repmat(obj.y_mem,size(x,1),1);
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        % train the algorithm
        function train(obj,~,y)
            obj.y_mem = y;
        end
    end
end
