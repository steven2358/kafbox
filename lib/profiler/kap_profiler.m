% Profiler extension for Kernel Affine Projection algorithm with Coherence
% Criterion
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef kap_profiler < kap
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
        prev_dict_size = 0; % previous dictionary size for growth check
    end
    
    methods

        function kaf = kap_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            kaf = kaf@kap(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if kaf.prev_dict_size < m, % growing
                m1 = m - 1;
                m2 = m - 1;
            else
                m1 = m;
                m2 = m;
            end
            m3 = size(kaf.memx,1);
            p2 = size(kaf.memx,1);
            
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
        
        % kernel(kaf.memx,kaf.dict,kaf.kerneltype,kaf.kernelpar);
        % kernel: m3*m2
        
        % kaf.alpha = kaf.alpha + ...
        %        kaf.eta*H'/...
        %        (kaf.eps*eye(size(H,1)) + H*H')*...
        %        (kaf.memd - H*kaf.alpha);
        % sum: m2 + p2^2*m2 - p2^2 + m2*(p2-1) + p2*m2 % without division
        % mult: m2*p2 + p2^2*m2 + p2*p2 +  p2*m2 % without division
        
        % matrix division in previous operation, assuming Gaussian elimination
        % sum: (p2-1)*p2*(2*p2-1)/6 - (p2-1)*p2/2
        % mult: (p2-1)*p2*(2*p2-1)/6 - (p2-1)*p2/2
        % div: p2 - 1
        
        %%
        
        function kaf = train_profiled(kaf,x,y)
            kaf.prev_dict_size = size(kaf.dict,1);
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m + m*size(kaf.dict,2) + kaf.p + kaf.p*size(kaf.memx,2)); % 8 bytes for double precision
            % alpha, dict, memd, memx
        end
        
    end
end
