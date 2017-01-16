% Demo: k-step ahead prediction on Lorenz attractor time-series data.
% This script trains and tests on separate splits of the data.

clear
close all;

%% PARAMETERS

k = 10; % prediction horizon
N_train = 1000; % number of training data
N_test = 1000; % number of test data, max 5000 - N_train
seed = 1; % random seed for reproducibility

% load data using helper function
[X,Y] = kafbox_data(struct('name','Lorenz','horizon',k,'embedding',6,'N',5000));

% Make a kernel adaptive filter object
kaf = krlst(struct('lambda',1,'M',100,'sn2',1E-4,'kerneltype','gauss','kernelpar',32));
% kaf = qklms(struct('eta',0.1,'epsu',.1,'kerneltype','gauss','kernelpar',32));

%% RUN ALGORITHM

% init random generator
rng('default')
rng(seed);

% corrupt outputs with noise
Y = Y + 3*randn(size(Y));

% split
X_train = X(1:N_train,:);
Y_train = Y(1:N_train);
X_test = X(N_train+1:N_train+N_test,:);
Y_test = Y(N_train+1:N_train+N_test);

Y_test_MSE = zeros(N_train,1);
for i=1:N_train,
    % train on data set 1
    if ~mod(i,floor(N_train/10)), fprintf('.'); end % progress indicator, 10 dots
    kaf.train(X(i,:),Y(i)); % train with one input-output pair
    
    % test on data set 2
    Y_test_est = kaf.evaluate(X_test);
    
    % store MSE
    SE = (Y_test-Y_test_est).^2; % out-of-sample error
    Y_test_MSE(i) = mean(SE);
end
fprintf('\n');

%% OUTPUT

% learning curve
figure
plot(10*log10(Y_test_MSE))
title('Learning curve')

% errors for final prediction
figure; hold all; plot(Y_test); plot(Y_test_est);
legend('original','prediction');
title(sprintf('%d-step ahead prediction result for %s on test set after %d training steps',...
    k,upper(class(kaf)),N_train));

fprintf('Final MSE: %.2f\n',Y_test_MSE(N_train));
