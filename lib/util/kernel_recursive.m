% Calculate a recursive kernel matrix.
% Author: Steven Van Vaerenbergh, 2016
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox

function K = kernel_recursive(X1,X2,ktype,kpar)

persistent k_prev;

switch ktype
        
    case 'gauss-recursive'
        sigma_i = kpar(1); % kernel width on input
        sigma = kpar(2); % kernel width on state
        
        if isempty(k_prev), % initialize persistent variable
            k_prev = 0;
        end
        
        K1 = kernel(X1,X2,'gauss',sigma_i);
        K = K1*exp((k_prev-1)/(2*sigma^2));
        
        k_prev = K(1); % store persistent variable
        
    case 'poly-recursive',
        p = kpar(1); % polynome order
        c = kpar(2); % additive constant
        sigma = kpar(3); % scaling for the state vector
        
        if isempty(k_prev), % initialize persistent variable
            k_prev = 0;
        end
        
        K = (X1*X2' + c + sigma*k_prev).^p;
        
        k_prev = K(1); % store persistent variable
        
    otherwise	% default case
        error ('unknown kernel type')
end
