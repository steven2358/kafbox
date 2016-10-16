% Demo: learn a sinc. Run one algorithm using its default parameters.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

close all
clear
rs = 1; randn('state',rs); rand('state',rs); %#ok<RAND>

%% PARAMETERS

N = 1000; % number of training data
N_test = 500; % number of test data
SNR = 20; % SNR in dB

algorithm = 'krls';

%% GENERATE DATA
x = randn(N,1);
x_test = linspace(min(x),max(x),N_test)';
y_ref = sinc([x;x_test]);
y = y_ref + sqrt(10^(-SNR/10)*var(y_ref))*randn(N+N_test,1);
y_test = y_ref(N+1:N+N_test);

%% RUN ALGORITHM
fprintf('%s: ',upper(algorithm));
Y_est = zeros(N_test,1);
kaf = feval(algorithm);
t1 = tic;
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end
    kaf.train(x(i),y(i));
end
y_est = kaf.evaluate(x_test);
MSE = mean((y_test-y_est).^2);

%% OUTPUT

fprintf(' %.2fs. MSE=%3.2fdB\n',toc(t1),10*log10(MSE))

figure; hold all
plot(x,y(1:N),'.')
    
plot(x_test,y_est,'LineWidth',2)
legend({'data',strrep(upper(algorithm),'_','-')})
axis([min(x)-0.5 max(x)+0.5 min(y)-0.5 max(y)+0.5]);
