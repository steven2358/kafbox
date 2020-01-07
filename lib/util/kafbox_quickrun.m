% Quick-run an algorithm on a dataset.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function kafbox_quickrun(kafname,datasetname,kafopt,dataopt)

rng(1)

t1 = tic;
if nargin > 2
    kaf = feval(kafname,kafopt);
else
    kaf = feval(kafname);
end

% get data
opt = struct('N_test',50);
[X,y,y_ref,X_test,y_test] = generate_channel_switch(opt); %#ok<ASGLU>

N = size(X,1);

N_switch = 500;

MSE = zeros(N,1);
for n=1:N
    if ~mod(n,floor(N/10))
        fprintf('.'); % progress indicator (10 dots)
    end
    
    y_est = kaf.evaluate(X_test);
    if n<=N_switch
        MSE(n) = mean((y_test(:,1)-y_est).^2);
    else
        MSE(n) = mean((y_test(:,2)-y_est).^2);
    end
    
    kaf.train(X(n,:),y(n)); % train with one input-output pair
end
fprintf(' %.2fs. Final MSE=%3.2fdB.\n',toc(t1),10*log10(mean(MSE(N-500:N,1))));

%% OUTPUT

figure;
plot(10*log10(MSE(:)))
legend(kafname)
