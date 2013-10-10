% Wrapper program for time series prediction. Performs progression of time
% and calls adaptive filter during each iteration.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function [e,kaf] = kafbox_predictionloop(kaf,X,y,vb)

if nargin<4,
    vb = 1;
end

N = size(X,1);
e = zeros(N,1);

for i=1:N,
    if ~mod(i,floor(N/10)) && vb,
        fprintf('.'); % progress indicator (10 dots)
    end
    
    y_est = kaf.evaluate(X(i,:)); % predict the next output
    e(i) = y(i)-y_est; % store error
    kaf = kaf.train(X(i,:),y(i)); % train with one input-output pair
end

if vb,
    fprintf('\n');
end;
