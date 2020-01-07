% Probabilistic Least-Mean-Squares Filter
%
% Jesus Fernandez-Bes, Victor Elvira, and Steven Van Vaerenbergh, "A
% probabilistic least-mean-squares filter," 2015 IEEE International
% Conference on Acoustics, Speech and Signal Processing (ICASSP), 2015.
% http://dx.doi.org/10.1109/ICASSP.2015.7178361
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef problms < linear_filter
    
    properties (GetAccess = 'public', SetAccess = 'private')
        sigma2_n = 1E-6; % variance of observation noise
        sigma2_d = 1E-4; % variance of filter weight diffusion
        lambda = 1; % forgetting factor
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        w = [];
        sigma2_k = 1E-2;
    end
    
    methods
        function obj = problms(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)'
                    if ismember(fn,fieldnames(obj))
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
            m = length(x);
            if numel(obj.w)==0 % initialize
                obj.w = zeros(m,1);
            end
            
            x_norm2 = x*x';
            K_k = (obj.sigma2_k + obj.sigma2_d) / ...
                ((obj.sigma2_k + obj.sigma2_d) * x_norm2 + obj.sigma2_n);
            
            err = y - obj.lambda*x*obj.w;
            obj.w = obj.lambda*obj.w + K_k*x'*err;
            obj.sigma2_k = (1 - K_k*x_norm2/m)*...
                (obj.sigma2_k + obj.sigma2_d);
        end
    end
end
