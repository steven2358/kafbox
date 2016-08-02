% Profiler extension for Recursive Least-Squares algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef rls_profiler < rls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
    end
    
    methods
        
        function obj = rls_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            obj = obj@rls(parameters);
        end
        
        function flops = lastflops(obj) % flops for last iteration
            m = size(obj.w,1);
            floptions = struct(...
                'sum', 2*m^2 - m + m + 2*m - 1 + 2*m^2 - m, ...
                'mult', m + m + m + 3*m^2, ...
                'div', 1 + 1);
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = obj.P*x'/(obj.lambda+x*obj.P*x');
        % sum: 2*m^2 - m
        % mult: m
        % div: 1
        
        % z = y - x*obj.w;
        % sum: m
        % mult: m
        
        % obj.w = obj.w + k*z;
        % sum: 2*m - 1
        % mult: m
        
        % obj.P = obj.lambda\(obj.P - k*x*obj.P);
        % sum: 2*m^2 - m
        % mult: 3*m^2
        % div: 1
        
        %%
        
        function obj = train_profiled(obj,x,y) % measures elapsed time of training
            t1 = tic;
            obj = obj.train(x,y);
            t2 = toc(t1);
            obj.elapsed = obj.elapsed + t2;
        end
        
        function bytes = lastbytes(obj) % bytes used in last iteration
            m = size(obj.w,1);
            bytes = 8*(m + m^2); % 8 bytes for double precision
            % w, P
        end
        
    end
end
