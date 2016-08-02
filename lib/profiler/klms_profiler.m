% Profiler extension for Kernel Least-Mean-Square algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef klms_profiler < klms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
        prev_mem_size = 0; % previous dictionary size for growth check
    end
    
    methods
        
        function kaf = klms_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            kaf = kaf@klms(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.mem,1);
            if kaf.prev_mem_size < m, % growing
                m1 = m - 1;
                floptions = struct(...
                    'sum', m1 - 1, ...
                    'mult', m1 + 1, ...
                    sprintf('%s_kernel',kaf.kerneltype), [m1, 1, size(kaf.mem,2)]);
                
                flops = kflops(floptions);
            else
                flops = 0;
            end
        end
        
        %% flops breakdown
        
        % k = kernel(kaf.mem,x,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % y_est = k'*kaf.alpha;
        % sum: m1 - 1
        % mult: m1
        
        % err = y - y_est;
        % sum: 1
        
        % kaf.alpha = [kaf.alpha; kaf.mu*err];
        % mult: 1
        
        %%
        
        function kaf = train_profiled(kaf,x,y)
            kaf.prev_mem_size = size(kaf.mem,1);
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.mem,1);
            bytes = 8*(m + m*size(kaf.mem,2)); % 8 bytes for double precision
            % alpha, mem
        end
        
    end
end
