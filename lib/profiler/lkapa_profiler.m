% Profiler extension for Leaky Kernel Affine Projection Algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef lkapa_profiler < lkapa
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
    end
    
    methods
        
        function kaf = lkapa_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            kaf = kaf@lkapa(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            % numbers of operations
            m1 = 1;
            m2 = 1;
            m3 = 1;
            m4 = 1;
            
            floptions = struct(...
                'sum', m1, ...
                'mult', m2, ...
                'div', m3, ...
                sprintf('%s_kernel',kaf.kerneltype), [m4,1,size(kaf.dict,2)]);
            
            flops = nan; % dummy. kflops(floptions);
        end
        
        %% flops breakdown
        
        % [space for comments on number of operations used above]
        
        %%
        
        function kaf = train_profiled(kaf,x,y)
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            m2 = 1;
            bytes = nan; % dummy. (m2 + m*size(kaf.dict,2)); % 8 bytes for double precision
            % [list variables]
        end
        
    end
end
