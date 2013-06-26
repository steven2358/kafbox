% Profiler extension for Kernel Least-Mean-Square algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef klms_profiler < klms
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            if kaf.grow
                m = size(kaf.dict,1);
                m1 = m - 1;
                floptions = struct(...
                    'sum', m1 - 1, ...
                    'mult', m1 + 1, ...
                    sprintf('%s_kernel',kaf.kerneltype), [m1, 1, size(kaf.dict,2)]);
                
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
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m + m*size(kaf.dict,2)); % 8 bytes for double precision
            % alpha, mem
        end
        
    end
end
