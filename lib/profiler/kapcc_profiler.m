% Profiler extension for Kernel Affine Projection algorithm with Coherence
% Criterion
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef kapcc_profiler < kapcc
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if kaf.growdict
                m1 = m - 1;
                m2 = m - 1;
            else
                m1 = m;
                m2 = m;
            end
            m3 = size(kaf.mem,1);
            p2 = size(kaf.mem,1);
            
            floptions = struct(...
                'sum', m2 + p2^2*m2 - p2^2 + m2*(p2-1) + p2*m2 + (p2-1)*p2*(2*p2-1)/6 - (p2-1)*p2/2, ...
                'mult', m2*p2 + p2^2*m2 + p2*p2 +  p2*m2 + (p2-1)*p2*(2*p2-1)/6 - (p2-1)*p2/2, ...
                'div', p2 - 1, ...
                sprintf('%s_kernel',kaf.kerneltype), [m1+m3*m2, 1, size(kaf.dict,2)]);
            
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % kernel(kaf.mem,kaf.dict,kaf.kerneltype,kaf.kernelpar);
        % kernel: m3*m2
        
        % kaf.alpha = kaf.alpha + ...
        %        kaf.eta*H'/...
        %        (kaf.eps*eye(size(H,1)) + H*H')*...
        %        (kaf.d - H*kaf.alpha);
        % sum: m2 + p2^2*m2 - p2^2 + m2*(p2-1) + p2*m2 % without division
        % mult: m2*p2 + p2^2*m2 + p2*p2 +  p2*m2 % without division
        
        % matrix division in previous operation, assuming Gaussian elimination
        % sum: (p2-1)*p2*(2*p2-1)/6 - (p2-1)*p2/2
        % mult: (p2-1)*p2*(2*p2-1)/6 - (p2-1)*p2/2
        % div: p2 - 1
        
        %%
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m + size(kaf.dict,2)) + 8*kaf.p(1 + size(kaf.mem,2)); % 8 bytes for double precision
        end
        
    end
end
