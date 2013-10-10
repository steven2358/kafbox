% Demo: learn a nonlinear channel with an abrupt change. Run and compare
% all algorithms using their default parameters.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

close all
clear all
% rs = 1; randn('state',rs); rand('state',rs); %#ok<RAND>

%% PARAMETERS

N = 1500; % number of training data
N_switch = 500; % iteration after which the channel switch occurs
N_test = 100; % number of test data
SNR = 30; % SNR in dB

%% PROGRAM

w = who;
for i=1:length(w), % copy all parameters to option structure
    eval(sprintf('opt.%s = %s;',w{i},w{i}))
end

% generate data
[X,y,y_ref,X_test,y_test] = generate_channel_switch(opt);

% get list of kernel adaptive filters in 'lib' folder
fdir = fileparts(which('kafbox_template.m'));
files = dir(fullfile(fdir,'*.m'));
[~,algorithms] = cellfun(@fileparts, {files.name}, 'UniformOutput',0);
for i=length(algorithms):-1:1,
    if ~exist(algorithms{i},'class'),
        algorithms(i) = []; % remove files that do not represent classes
    end
end

% perform online learning for each algorithm
fprintf('\n')
num_alg = length(algorithms);
titles = cell(num_alg,1);
MSE = zeros(N,num_alg);
MSE_final = zeros(1,num_alg);
for algo_ind=1:num_alg
    t1 = tic;
    algorithm = algorithms{algo_ind};
    fprintf('%2d. %9s: ',algo_ind,upper(algorithm));
    titles{algo_ind} = strrep(upper(algorithm),'_','\_');

    kaf = feval(algorithm);
    for i=1:N,
        if ~mod(i,floor(N/10)), fprintf('.'); end
        
        y_est = kaf.evaluate(X_test);
		if i<=N_switch
			MSE(i,algo_ind) = mean((y_test(:,1)-y_est).^2);
		else
			MSE(i,algo_ind) = mean((y_test(:,2)-y_est).^2);
		end			
        
        kaf = kaf.train(X(i,:),y(i));
    end
    MSE_final(algo_ind) = mean(MSE(N-500:N,algo_ind));
    
    fprintf(' %.2fs. Final MSE=%3.2fdB\n',toc(t1),10*log10(MSE_final(algo_ind)))
end

%% OUTPUT

% plot results in different "leagues"
[MSE_final_sorted,ind] = sort(MSE_final,'descend');
num_fig = ceil(num_alg/5);
counter = [rem(num_alg,5) rem(num_alg,5)+1];
for fig_ind=num_fig:-1:1
    figure; hold all
    for i=1:counter(1),
        counter(2) = counter(2) - 1;
        plot(10*log10(MSE(:,ind(counter(2)))),'LineWidth',1)
    end
    title(sprintf('League %d',fig_ind))
    legend(titles(ind(counter(2)+counter(1)-1:-1:counter(2))))
    counter(2) = counter(2)+counter(1)+5;
    counter(1) = 5;
    axis([0 N 5*floor(min(10*log10(MSE(:)))/5) 0]);
end
