% Profiler extension for Kernel Normalized Least-Mean-Square algorithm with
% Coherence Criterion
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef knlms_profiler < knlms
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if kaf.grow
                m1 = m - 1;
                n2 = 1;
            else
                m1 = m;
                n2 = 0;
            end
            m3 = m;
            floptions = struct(...
                'sum', 3*m3, ...
                'mult', 3*m3 + 1, ...
                'kernel', [kaf.kerneltype, m1+n2, size(kaf.dict,2)]);
            
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % h = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
        % kernel: n2 % 1 element when growing, 0 when not
        
        % kaf.alpha = kaf.alpha + kaf.eta / (kaf.eps + h*h') * (y - h*kaf.alpha) * h';
        % sum: 3*m3
        % mult: 3*m3 + 1
        % div: 1
        
        %%
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m + m*size(kaf.dict,2)); % 8 bytes for double precision
        end
        
    end
end
