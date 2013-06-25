% Profiler extension for Extended Kernel Recursive Least Squares algorithm
% Linear Dependency
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef exkrls_profiler < exkrls
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.mem,1);
            if kaf.grow
                m1 = m;
                m2 = m - 1;
                floptions = struct(...
                    'sum', 1 + m2 - 1 + m2 + 1 + m2 + 2*m2 - 1 + 1 + m2^2 + m2 - 1, ...
                    'mult', m2 + m2 + 2 + m2 + m2 + 2 + 4 + 2*m2^2 + m2 + 2, ...
                    'div', 1 + 1, ...
                    sprintf('%s_kernel',kaf.kerneltype), [m1, 1, size(kaf.dict,2)]);
            else
                floptions = struct('sum', 1); % no kernel calculation
            end
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = kernel([kaf.mem; x],x,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % kaf.i = kaf.i + 1;
        % sum: 1
        
        % z = kaf.Q*kt;
        % sum: m2 - 1
        % mult: m2
        
        % r = kaf.beta^kaf.i*kaf.rho + ktt - kt'*z;
        % sum: m2 + 1
        % mult: m2 + 2 % assumes previous kaf.beta^kaf.i is stored
        
        % err = y - kt'*kaf.alpha;
        % sum: m2
        % mult: m2
        
        % kaf.alpha = kaf.alphaf*[kaf.alpha - z*err/r; err/r];
        % sum: 2*m2 - 1
        % mult: m2 + 2
        % div: 1
        
        % dummy = kaf.alphaf^2 + kaf.beta^kaf.i*kaf.q*kaf.rho;
        % sum: 1
        % mult: 4
        
        % kaf.rho = kaf.rho/dummy;
        % div: 1
        
        % kaf.Q = kaf.alphaf^2/(r*dummy)*[kaf.Q*r + z*z', -z; -z', 1];
        % sum: m2^2 + m2 - 1
        % mult: 2*m2^2 + m2 + 2
        % div: 1
        
        %%
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m^2 + m + 1 + m*size(kaf.dict,2)); % 8 bytes for double precision
        end
        
    end
end
