% Base class for all estimators in KAFBOX.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

##classdef base_estimator < matlab.mixin.Copyable
classdef base_estimator < handle
    methods
        % get parameter names for the estimator
        function names = get_param_names(obj)
            names = fieldnames(obj);
        end

        % get parameters
        function params = get_params(obj)
            params = struct;
            for fn = fieldnames(obj)'
                params.(fn{1}) = obj.(fn{1});
            end
        end

        % set parameters
        function set_params(obj,params)
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(params)'
                    if ismember(fn{1},fieldnames(obj))
                        values = params.(fn{1});
                        obj.(fn{1}) = values;
                    else
                        warning('Unknown parameter: %s.',fn{1});
                    end
                end
            end
        end
    end
end
