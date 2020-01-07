% Reproduces an experiment similar to the one from figure 2 in "A
% sliding-window kernel RLS algorithm and its application to nonlinear
% channel identification".
%
% MSE performance of sliding-window kernel RLS for two different window
% sizes. Execution time: < 1 minute (Intel Pentium Core2 Duo).
%
% S. Van Vaerenbergh, J. Via, and I. Santamaria, "A Sliding-Window Kernel
% RLS Algorithm and Its Application to Nonlinear Channel Identification,"
% 2006 IEEE International Conference on Acoustics, Speech and Signal
% Processing (ICASSP), May 2006,
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab
% https://github.com/steven2358/kafbox/

clear
close all

%% PARAMETERS
opt.chlen = 4; % channel length, and input embedding
opt.N = 1500; % number of training data
opt.N_switch = 500; % iteration after which the channel switch occurs
opt.N_test = 100; % number of test data
opt.SNR = 20; % SNR in dB

num_sim = 10;

setups{1} = struct('M',75,'c',1E-2,'kerneltype','gauss','kernelpar',5);
setups{2} = struct('M',150,'c',1E-2,'kerneltype','gauss','kernelpar',5);

titles = {'SW-KRLS, M=75','SW-KRLS, M=150'};

%% PREPARE DATA
fprintf('Fig. 2 from "A sliding-window kernel RLS algorithm and its \n');
fprintf('application to nonlinear channel identification".\n');

% perform online learning for each algorithm
fprintf('\n')
N = opt.N;
N_switch = opt.N_switch;
num_alg = length(setups);
MSE = zeros(N,num_alg);
MSE_final = zeros(1,num_alg);
for sim_ind = 1:num_sim
    % generate data
    [X,y,y_ref,X_test,y_test] = generate_channel_switch(opt);
    
    fprintf('SIM %d:\n',sim_ind)
    for algo_ind = 1:num_alg
        t1 = tic;
        
        kaf = swkrls(setups{algo_ind});
        
        fprintf('%9s M=%3d: ',upper(class(kaf)),kaf.M);
        
        mse = zeros(N,1);
        for i=1:N
            if ~mod(i,floor(N/10)), fprintf('.'); end
            
            y_est = kaf.evaluate(X_test);
            if i<=N_switch
                mse(i) = mean((y_test(:,1)-y_est).^2);
            else
                mse(i) = mean((y_test(:,2)-y_est).^2);
            end
            
            kaf.train(X(i,:),y(i));
        end
        MSE(:,algo_ind) = MSE(:,algo_ind) + mse/num_sim;
        mse_final = mean(mse(N-500:N));
        
        fprintf(' %.2fs. Final MSE=%3.2fdB\n',toc(t1),...
            10*log10(mse_final))
    end
    fprintf('\n');
end
fprintf('\n');

%% OUTPUT

figure; hold all
for algo_ind=1:num_alg
    plot(10*log10(MSE(:,algo_ind)),'LineWidth',1)
    
    axis([0 N 5*floor(min(10*log10(MSE(:)))/5) 0]);
end
xlabel('iteration')
ylabel('MSE (dB)')
legend(titles)
    