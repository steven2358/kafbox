% Profiler extension for Kernel Least-Mean-Square algorithm
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef klms_profiler < klms
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            m1 = m - 1; % only valid when growing, otherwise small error of 1 flop in summation only
            floptions = struct(...
                'sum', m1 - 1, ...
                'mult', m1 + 1, ...
                'kernel', [kaf.kerneltype,m1,size(kaf.dict,2)]);
            
            flops = kflops(floptions);
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
        end
        
    end
end
