% Partly reproduces figure 3 from "Kernel Recursive Least-Squares Tracker
% for Time-Varying Regression". (Only 3 algorithms, only 1 MC simulation.)
%
% MSE performance comparison of different tracking algorithms on a
% communications channel that shows an abrupt change at iteration 500.
% Execution time: < 1 minute (Pentium Core2 Duo).
%
% S. Van Vaerenbergh, M. Lazaro-Gredilla, and I. Santamaria, "Kernel
% Recursive Least-Squares Tracker for Time-Varying Regression," IEEE
% Transactions on Neural Networks and Learning Systems, vol. 23, no. 8, pp.
% 1313-1326, Aug. 2012, http://dx.doi.org/10.1109/TNNLS.2012.2200500
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab
% http://sourceforge.net/projects/kafbox/

clear all
close all

%% PARAMETERS

N = 1500; % number of training data points
N_test = 500; % number of test data points
Nswitch = 500; % switch from model 1 to model 2 after Nswitch iterations
B1 = [1.0000, -0.3817, -0.1411, 0.5789, 0.191]; % model 1 linear filter
B2 = [1.0000, -0.0870, 0.9852, -0.2826, -0.1711]; % model 2 linear filter
f = inline('tanh(x)'); % Wiener system nonlinearity
SNR = 20; % SNR in dB

embedding = 5; % time-embedding

setups{1} = swkrls(struct('c',1E-2,'M',50,'kerneltype','gauss','kernelpar',1));
setups{2} = krlst(struct('lambda',.999,'M',50,'sn2',1E-2,'kerneltype','gauss','kernelpar',1));
setups{3} = norma(struct('lambda',1E-2,'tau',1500,'eta',0.1,'kerneltype','gauss','kernelpar',1)); % large tau -> very slow

numsim = 1;

%% PREPARE DATA

fprintf('Fig. 3 from "Kernel Recursive Least-Squares Tracker for\n');
fprintf('Time-Varying Regression".\n')

% generate Gaussian input data s
s = rand(N+N_test,1);
s_mem = zeros(N+N_test,embedding);
for i = 1:embedding,
    s_mem(i:N+N_test,i) = s(1:N+N_test-i+1);	% time-embedding
end
s = s_mem(1:N+N_test,:);	% input data, stored in columns
s_train = s_mem(1:N,:);	% input train data, stored in columns
s_test = s_mem(N+1:N+N_test,:);	% input test data, stored in columns

% generate internal data x and output y
X1 = s_mem(1:Nswitch,:)*B1';
X2 = s_mem(Nswitch+1:N,:)*B2';
X = [X1;X2];
Y_nn = f(X); % noiseless Y
vary = var(Y_nn);
noisevar = 10^(-SNR/10)*vary;
noise = sqrt(noisevar)*randn(N,1);
Y = Y_nn + noise;		% noisy output data

X_test1 = s_mem(N+1:N+N_test,:)*B1';
X_test2 = s_mem(N+1:N+N_test,:)*B2';
noise_test1 = sqrt(noisevar)*randn(N_test,1);
noise_test2 = sqrt(noisevar)*randn(N_test,1);
Y_test1 = f(X_test1) + noise_test1;	% noisy output test data, model 1
Y_test2 = f(X_test2) + noise_test2;	% noisy output test data, model 2

%% RUN ALGORITHMS
t1 = tic;

num_setup = length(setups);
MSE = zeros(N,num_setup);
titles = cell(num_setup,1);

for setup_ind=1:length(setups)
    t2 = tic;
    kaf = setups{setup_ind};
    
    titles{setup_ind} = upper(class(kaf));
    fprintf('%s\t',titles{setup_ind});
    for n=1:N,
        if ~mod(n,round(N/10)), fprintf('.'); end
        if n<=Nswitch,
            Y_test = Y_test1;
        else
            Y_test = Y_test2;
        end
        Y_est = kaf.evaluate(s_test); % test on test set
        err = Y_test - Y_est;
        MSE(n,setup_ind) = mean(err.^2);
        
        kaf = kaf.train(s_train(n,:),Y(n)); % train with one input-output pair
    end
    fprintf(' %.2f seconds\n',toc(t2));
end

toc(t1)
%% OUTPUT

figure
plot(10*log10(MSE),'LineWidth',2); grid on
xlabel('iteration')
ylabel('MSE (dB)')
legend(titles)
