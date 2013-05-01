% Profiler extension for Sliding-Window Kernel Recursive Least Squares.
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef swkrls_profiler < swkrls
    
    methods
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if ~kaf.prune, m = m-1; end
            floptions = struct(...
                'sum',3*m^2-m+1 + kaf.prune*(m^2),...
                'mult',3*m^2+m+1 + kaf.prune*(m^2+m),...
                'div',1 + kaf.prune*(1),...
                'kernel',[kaf.kerneltype,m,size(x,2)]);
            flops = kflops(floptions);
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*m*(m+2+size(x,2)); % 8 bytes for double precision
        end
    end
    
end
