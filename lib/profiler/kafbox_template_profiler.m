% Profiler extension for Kernel Recursive Least-Squares Tracker algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef kafbox_template_profiler < kafbox_template
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
    end
    
    methods
        
        function kaf = kafbox_template_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            kaf = kaf@kafbox_template(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            % numbers of operations
            m1 = 1;
            m2 = 1;
            m3 = 1;
            m4 = 1;
            
            floptions = struct(...
                'sum', m1, ...
                'mult', m2, ...
                'div', m3, ...
                sprintf('%s_kernel',kaf.kerneltype), [m4,1,size(kaf.dict,2)]);
            
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % space for comments on number of operations used above
        
        %%
        
        function kaf = train_profiled(kaf,x,y)
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m^2 + m^2 + m + 2 + m*size(kaf.dict,2)); % 8 bytes for double precision
            % Q, Sigma, mu, nums02ML, dens02ML, dict
        end
        
    end
end
