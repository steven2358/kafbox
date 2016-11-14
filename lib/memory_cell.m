% Dummy prediction method that always returns the last seen output.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef memory_cell < handle
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        y_mem;
    end
    
    methods
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if ~isempty(kaf.y_mem)
                y_est = repmat(kaf.y_mem,size(x,1),1);
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function train(kaf,~,y) % train the algorithm
            kaf.y_mem = y;
        end
    end
end
