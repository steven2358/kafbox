% Demo: k-step ahead prediction on Lorenz attractor time-series data

clear
close all;

%% PARAMETERS

k = 10; % prediction horizon

% load data using helper function
[X,Y] = kafbox_data(struct('name','Lorenz','horizon',k,...
    'embedding',6,'N',5000));

% Make a kernel adaptive filter object
kaf = krlst(struct('lambda',1,'M',100,'sn2',1E-6,'kerneltype','gauss','kernelpar',32));
% kaf = qklms(struct('eta',0.5,'epsu',.1,'kerneltype','gauss','kernelpar',32));

%% RUN ALGORITHM
N = size(X,1);
Y_est = zeros(N,1);
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end % progress indicator, 10 dots
    Y_est(i) = kaf.evaluate(X(i,:)); % predict the next output
    kaf.train(X(i,:),Y(i)); % train with one input-output pair
end
fprintf('\n');
SE = (Y-Y_est).^2; % test error

%% OUTPUT
fprintf('MSE after first 1000 samples: %.2fdB\n\n',...
    10*log10(mean(SE(1001:end))));

figure; hold all; plot(Y); plot(Y_est);
legend('original','prediction');
title(sprintf('%d-step ahead prediction %s on Lorenz time series',...
    k,upper(class(kaf))));
