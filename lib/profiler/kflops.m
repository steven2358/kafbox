function fl = kflops(operations)
% Calculate the number of FLOPS needed to perform the specified operations.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

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
        case 'gauss_kernel'
            n1 = operations.(op{1})(1);
            n2 = operations.(op{1})(2);
            m = operations.(op{1})(3);
            
            fl = fl + kflops(struct('sum',n1*n2*(2*m-1),...
                'mult',n1*n2*(m+1),'exp',n1*n2));
            
            % breakdown:
            % fl = kaf_flops('sum',N1*N2*M) + ... % all X1(i,:)-X(2(j,:)
            % kaf_flops('mult',N1*N2*M) + ... % square all elements
            % kaf_flops('sum',N1*N2*(M-1)) + ... % sum M-1 elements per entry
            % kaf_flops('mult',N1*N2) + ... % multiply by 1/(2*sgm^2)
            % kaf_flops('exp',N1*N2);
            
        otherwise
            error('unknown option')
    end
end
