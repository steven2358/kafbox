% Profiler extension for Naive Online regularized Risk Minimization
% Algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef norma_profiler < norma
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
        prev_mem_size = 0; % previous memory size for prune check
    end
    
    methods
        
        function kaf = norma_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            kaf = kaf@norma(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            
            m = size(kaf.mem,1);
            if kaf.prev_mem_size < m, % growing (no pruning)
                m1 = m-1;
            else
                m1 = m;
            end
            floptions = struct(...
                'sum', m1 - 1 + 1 + 1, ...
                'mult', 2*m1 + 3 + 1, ...
                'exp', 1, ...
                sprintf('%s_kernel',kaf.kerneltype), [m1,1,size(kaf.mem,2)]);
            
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = kernel(kaf.mem,x,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % y_est = k'*(kaf.alpha.*kaf.beta(length(kaf.alpha):-1:1));
        % sum: m1 - 1
        % mult: 2*m1
        
        % err = y - y_est;
        % sum: 1
        
        % kaf.alpha = (1-kaf.lambda*kaf.eta*kaf.t^kaf.tcoeff)*kaf.alpha;
        % sum: 1
        % mult: 3
        % exp: 1
        
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
            bytes = 4 + 8*(m + m*size(kaf.mem,2)); % 8 bytes for double precision, 4 for uint32
            % t, alpha, mem
        end
        
    end
end
