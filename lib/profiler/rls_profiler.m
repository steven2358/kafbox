% Profiler extension for Recursive Least-Squares algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef rls_profiler < rls
    
    properties (GetAccess = 'public', SetAccess = 'public')
        elapsed = 0;
    end
    
    methods
        
        function kaf = rls_profiler(parameters) % constructor
            kaf = kaf@rls(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.w,1);
            floptions = struct(...
                'sum', 2*m^2 - m + m + 2*m - 1 + 2*m^2 - m, ...
                'mult', m + m + m + 3*m^2, ...
                'div', 1 + 1);
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = kaf.P*x'/(kaf.lambda+x*kaf.P*x');
        % sum: 2*m^2 - m
        % mult: m
        % div: 1
        
        % z = y - x*kaf.w;
        % sum: m
        % mult: m
        
        % kaf.w = kaf.w + k*z;
        % sum: 2*m - 1
        % mult: m
        
        % kaf.P = kaf.lambda\(kaf.P - k*x*kaf.P);
        % sum: 2*m^2 - m
        % mult: 3*m^2
        % div: 1
        
        %%
        
        function kaf = train_elapsed(kaf,x,y) % measures elapsed time of training
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.w,1);
            bytes = 8*(m + m^2); % 8 bytes for double precision
        end
        
    end
end
