% Almost reproduces figure 1a from "Multikernel Adaptive Filtering".
%
% Comparison of the performance of LMS, KNLMS and MKNLMS-CS in nonlinear 
% channel equalization. Execution time: < 1 minute (Intel Pentium Core2
% Duo).
% 
% Masahiro Yukawa, "Multikernel Adaptive Filtering", IEEE Transactions on
% Signal Processing, vol.60, no.9, pp.4672-4682, Sept. 2012.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab
% https://github.com/steven2358/kafbox/

clear
close all

%% PARAMETERS
N_tr = 3000; % number of training data
N_te = 500; % number of testing data
SNR = 20;

embedding = 5; % input embedding
delay = 2; % equalization delay

mk_thresh = 0.36; % threshold value from paper
% mk_thresh = 0.9; % alternative value

mk_kernelpars = [.2 2]; % kernel parameter values from paper
% mk_kernelpars = [.2 .5 1 2 5 10 20 50]; % alternative values

% remark 1 from paper determines mu0
setups{1} = lms(struct('mu',.01));
setups{2} = knlms(struct('mu0',mk_thresh^(.2/min(mk_kernelpars)),'eta',0.2,'eps',3E-2,'kerneltype','gauss','kernelpar',1/sqrt(.2)));
setups{3} = knlms(struct('mu0',mk_thresh^(.5/min(mk_kernelpars)),'eta',0.2,'eps',3E-2,'kerneltype','gauss','kernelpar',1/sqrt(.5)));
setups{4} = mknlms_cs(struct('delta',mk_thresh,'eta',0.2,'rho',6E-2,'kerneltype','gauss','kernelpars',1./sqrt(mk_kernelpars)));

%% PREPARE DATA

fprintf('Fig. 1a from "Multikernel Adaptive Filtering".\n')

u = randn(N_tr+N_te+embedding-1,1)>0;
u = 2*u-1; % binary input

z = u + 0.5*[0;u(1:end-1)]; % output of linear channel
varz = var(z);
noisevar = 10^(-SNR/10)*varz;
ns = sqrt(noisevar)*randn(length(u),1); % channel noise
y = z - 0.9*z.^2 + ns; % output of nonlinear channel

X_all = zeros(N_tr+N_te,embedding); % time-embedding
for k=1:embedding,
    X_all(:,k) = y(k:N_tr+N_te+k-1);
end

X_tr = X_all(1:N_tr,:); % training input data
X_te = X_all(N_tr+1:N_tr+N_te,:); % test input data

T_tr = u(delay:delay+N_tr-1); % training desired output
T_te = u(delay+N_tr:delay+N_tr+N_te-1); % training desired output

%% RUN ALGORITHMS
tic

num_setup = length(setups);
MSE = zeros(N_tr,num_setup);

for setup_ind = 1:num_setup,
    kaf = setups{setup_ind};
    
    for n=1:N_tr
        if ~mod(n,floor(N_tr/10)), fprintf('.'); end % progress indicator, 10 dots
        
        t_te = kaf.evaluate(X_te); % test on test set
        err = T_te - t_te;
        MSE(n,setup_ind) = mean(err.^2);
        
        kaf = kaf.train(X_tr(n,:),T_tr(n)); % train with one input-output pair
    end
    if setup_ind == 1
        fprintf('\n');
    else
        fprintf(' Dict. size: %d\n',size(kaf.dict,1));
    end
end

toc
%% OUTPUT

figure
semilogy(MSE,'LineWidth',2); grid on

legend('LMS','KNLMS 1','KNLMS 2','MKNLMS-CS')
xlabel('Iteration number')
ylabel('MSE')

