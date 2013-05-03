function fl = kflops(operations)
% Calculate the number of FLOPS needed to perform the specified operations.
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

% values for x86 processor, change as desired
flops_sum = 1;
flops_mult = 1;
flops_div = 8;
flops_sqrt = 8;
flops_exp = 20;

fl = 0;
f = fields(operations);
for op = f'
    num = operations.(op{1});
    switch op{1}
        case 'sum'
            fl = fl + num*flops_sum;
        case 'mult'
            fl = fl + num*flops_mult;
        case 'div'
            fl = fl + num*flops_div;
        case 'exp'
            fl = fl + num*flops_exp;
        case 'sqrt'
            fl = fl + num*flops_sqrt;
    end
end
