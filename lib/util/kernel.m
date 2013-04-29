% Calculate the kernel matrix for two data sets.
% Author: Steven Van Vaerenbergh, 2013
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function K = kernel(X1,X2,ktype,kpar)
N1 = size(X1,1);
N2 = size(X2,1);

switch ktype
    case 'gauss' % RBF kernel
        norms1 = sum(X1.^2,2);
        norms2 = sum(X2.^2,2);
        
        mat1 = repmat(norms1,1,N2);
        mat2 = repmat(norms2',N1,1);
        
        dist2 = mat1 + mat2 - 2*X1*X2';	% full distance matrix
        K = exp(-dist2/(2*kpar^2));
        
    case 'gauss-diag' % diagonal of RBF kernel
        K = exp(-sum((X1-X2).^2,2)/(2*kpar^2));
        
    case 'poly'	% polynomial kernel
        p = kpar(1); % polynome order
        c = kpar(2); % additive constant
        
        K = (X1*X2' + c).^p;
        
    case 'linear' % linear kernel
        K = X1*X2';
        
    otherwise	% default case
        error ('unknown kernel type')
end
