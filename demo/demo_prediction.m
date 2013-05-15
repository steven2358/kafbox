%% PARAMETERS
datafile = 'lorenz.dat'; % Lorenz attractor time-series data
L = 6; % length of time-embedding
horizon = 1; % 1-step ahead prediction

% make a kernel adaptive filter object of class aldkrls with options: 
% ALD threshold 1E-4, Gaussian kernel, and kernel width 32
kaf = aldkrls(struct('nu',1E-4,'kerneltype','gauss','kernelpar',32));

%% PREPARE DATA
data = load(datafile); % 1-dimensional time-series
data = data(:);
N = length(data)-horizon; % number of data to use (all data)
X = zeros(N,L);
for i = 1:L,
    X(i:N,i) = data(1:N-i+1); % time embedding of length L
end
Y = data(1+horizon:N+horizon); % vector with desired outputs

%% RUN ALGORITHM
Y_est = zeros(N,1);
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end % progress indicator, 10 dots
    Y_est(i) = kaf.evaluate(X(i,:)); % predict the next output
    kaf = kaf.train(X(i,:),Y(i)); % train with one input-output pair
end
fprintf('\n');
SE = (Y-Y_est).^2; % test error

%% OUTPUT
fprintf('MSE after first 1000 samples: %.2fdB\n\n',10*log10(mean(SE(1001:end))));
