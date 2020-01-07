% 1-step ahead prediction on Mackey Glass time series.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

close all; clear

%% PARAMETERS

h = 1; % prediction horizon
L = 10; % embedding
n_train = 200; % train data
n_test = 100; % test data

sigma = 1; % kernel parameter
lms_lr = 0.2; % lms learning rate
klms_lr = 0.2; % klms learning rate
krls_nu = 1E-4; % krls precision parameter

% select algorithms
i=1;
% algos{i} = lms(struct('mu',lms_lr)); i=i+1;
algos{i} = klms(struct('M',Inf,'kerneltype','gauss','kernelpar',sigma,'eta',klms_lr)); i=i+1;
% algos{i} = krls(struct('nu',krls_nu,'kerneltype','gauss','kernelpar',sigma)); i=i+1;

%% PREPARE DATA

% load data
[X,y] = kafbox_data(struct('name','mg30','embedding',L,'horizon',h,...
    'N',n_train+n_test));

% split data
X_train = X(1:n_train,:);
y_train = y(1:n_train);
X_test = X(n_train+1:n_train+n_test,:);
y_test = y(n_train+1:n_train+n_test);

%% RUN ALGORITHMS
n_algos = length(algos);
MSE = zeros(n_train,n_algos);
y_est_final = zeros(n_test,n_algos);
titles = cell(n_algos,1);
for j=1:n_algos
    kaf = algos{j};
    titles{j} = upper(class(kaf)); % store algorithm name
    fprintf('Training %s',titles{j})
    for i=1:n_train
        if ~mod(i,floor(n_train/10)), fprintf('.'); end % progress indicator
        
        kaf.train(X_train(i,:),y_train(i)); % train with one input-output pair
        
        y_est = kaf.evaluate(X_test); % evaluate on test data
        MSE(i,j) = mean((y_test-y_est).^2); % test error
    end
    fprintf('\n');
    y_est_final(:,j) = y_est;
end

%% OUTPUT
figure;
plot(10*log10(MSE));
title('Learning curves')
legend(titles);
xlabel('iteration'); ylabel('MSE (db)');

fprintf('Final MSE: %.2fdB\n\n',10*log10(MSE(end)));

figure; hold all;
plot(y_test,'LineWidth',2);
titles2 = {'original'};
line_styles = {'--','-.'};
for j=1:n_algos
    plot(y_est_final(:,j),'LineWidth',2,'LineStyle',line_styles{j});
end
titles2(2:n_algos+1) = titles;
legend(titles2);
title(sprintf('%d-step ahead prediction on Mackey Glass time series',h));

% export_fig('test.pdf')
