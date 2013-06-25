% Profiler extension for Kernel Least-Mean-Square algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef qklms_profiler < qklms
    
    methods
        
        function kaf = qklms_profiler(parameters) % constructor
            kaf = kaf@qklms(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.mem,1);
            
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
                m2 = m1; m3 = size(kaf.mem,2);
            end
            
            floptions = struct(...
                'sum', m1 - 1 + 1 + (2*m2-1)*m3 + n4 + n5, ...
                'mult', m1 + m2*m3 + n4 + n5, ...
                sprintf('%s_kernel',kaf.kerneltype), [m1,1,size(kaf.mem,2)]);
            
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
        
        function kaf = train_elapsed(kaf,x,y) % measures elapsed time of training
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.mem,1);
            bytes = 8*(m + m*size(kaf.mem,2)); % 8 bytes for double precision
        end
        
    end
end
