% Profiler extension for Naive Online regularized Risk Minimization
% Algorithm
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef norma_profiler < norma
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if ~kaf.prune
                m1 = m-1;
            else
                m1 = m;
            end
            floptions = struct(...
                'sum', m1 - 1 + 1, ...
                'mult', 2*m1 + 1, ...
                'kernel', [kaf.kerneltype,m1,size(kaf.dict,2)]);
            
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
        
        % kaf.alpha = [kaf.alpha; kaf.mu*err];
        % mult: 1
        
        %%
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m + m*size(kaf.dict,2) + kaf.tau); % 8 bytes for double precision
        end
        
    end
end
