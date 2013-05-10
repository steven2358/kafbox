% Profiler extension for Kernel Least-Mean-Square algorithm
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef qklms_profiler < qklms
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            
            if kaf.grow
                m1 = m-1;
                n4 = 1;
                n5 = 1;
            else
                m1 = m;
                n4 = 1;
                n5 = 0;
            end
            
            if strcmp('kaf.kernel','gauss')
                m2 = 0; m3 = 0; % quantization criterion can use values calculated by kernel
            else
                m2 = m1; m3 = size(kaf.dict,2);
            end
            
            floptions = struct(...
                'sum', m1 - 1 + 1 + (2*m2-1)*m3 + n4 + n5, ...
                'mult', m1 + m2*m3 + n4 + n5, ...
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
        
        % [d2,j] = min(sum((kaf.mem - repmat(x,m,1)).^2,2));
        % sum: (2*m2-1)*m3
        % mult: m2*m3
        
        % kaf.alpha(j) = kaf.alpha(j) + kaf.eta*err;
        % sum: n4
        % mult: n5
        
        % kaf.alpha = [kaf.alpha; kaf.mu*err];
        % mult: n4
        
        %%
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m + m*size(kaf.dict,2)); % 8 bytes for double precision
        end
        
    end
end
