% Profiler extension for Kernel Recursive Least Squares with Approximate
% Linear Dependency
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef aldkrls_profiler < aldkrls
    
    methods
        
        function kaf = aldkrls_profiler(parameters) % constructor
            kaf = kaf@aldkrls(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if kaf.grow
                m1 = m;
                m2 = m - 1;
                m3 = m - 1;
                floptions = struct(...
                    'sum', m2^2 - m2 + m2 + m3^2 + m3 + m3, ...
                    'mult', m2^2 + m2 + m3^2 + m3 + m3 + 1 + m3, ...
                    'div', 1, ...
                    sprintf('%s_kernel',kaf.kerneltype), [m1,1,size(kaf.dict,2)]);
            else
                m1 = m + 1;
                m2 = m;
                m4 = m;
                floptions = struct(...
                    'sum', m2^2 - m2 + m2 + m4^2 + 2*m4^2 - m4 + m4^2 + m4, ...
                    'mult', m2^2 + m2 + 2*m4^2 + m4 + 2*m4^2 + m4^2 + 2*m4, ...
                    'div', 1, ...
                    sprintf('%s_kernel',kaf.kerneltype), [m1,1,size(kaf.dict,2)]);
            end
            
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % k = kernel([kaf.dict; x],x,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % at = kaf.Kinv*kt; % check linear dependency
        % sum: m2^2 - m2
        % mult: m2^2
        
        % delta = ktt - kt'*at; % check linear dependency
        % sum: m2
        % mult: m2
        
        % kaf.Kinv = 1/delta*[delta*kaf.Kinv + at*at', -at; -at', 1]; % grow Kinv
        % sum: m3^2
        % mult: m3^2 + m3
        % div: 1
        
        % ode = 1/delta*(y-kt'*kaf.alpha); % grow Kinv
        % sum: m3
        % mult: m3 + 1
        
        % kaf.alpha = [kaf.alpha - at*ode; ode]; % grow Kinv
        % sum: m3
        % mult: m3
        
        % q = kaf.P*at/(1+at'*kaf.P*at); % only update alpha
        % sum: m4^2
        % mult: 2*m4^2 + m4
        % div: 1
        
        % kaf.P = kaf.P - q*(at'*kaf.P); % only update alpha
        % sum: 2*m4^2 - m4
        % mult: 2*m4^2
        
        % kaf.alpha = kaf.alpha + kaf.Kinv*q*(y-kt'*kaf.alpha); % only update alpha
        % sum: m4^2 + m4
        % mult: m4^2 + 2*m4
        
        %%
        
        function kaf = train_elapsed(kaf,x,y) % measures elapsed time of training
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*m*(2*m + 1 + size(kaf.dict,2)); % 8 bytes for double precision
        end
        
    end
end
