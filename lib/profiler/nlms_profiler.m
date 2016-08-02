% Profiler extension for Normalized Least Mean Squares algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef nlms_profiler < nlms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
    end
    
    methods
        
        function obj = nlms_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            obj = obj@nlms(parameters);
        end
        
        function flops = lastflops(obj) % flops for last iteration
            m = size(obj.w,1);
            floptions = struct(...
                'sum', 3*m, ...
                'mult', 3*m, ...
                'div', 1);
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % obj.w = obj.w + obj.mu/(obj.eps + x*x')*x*(y-x*obj.w');
        % sum: 3*m
        % mult: 3*m
        % div: 1
        
        %%
        
        function obj = train_profiled(obj,x,y)
            t1 = tic;
            obj = obj.train(x,y);
            t2 = toc(t1);
            obj.elapsed = obj.elapsed + t2;
        end
        
        function bytes = lastbytes(obj) % bytes used in last iteration
            m = size(obj.w,1);
            bytes = 8*m; % 8 bytes for double precision
            % w
        end
        
    end
end
