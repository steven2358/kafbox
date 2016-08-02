% Reproduces figure 2.12 from "Kernel Adaptive Filtering: A Comprehensive
% Introduction".
%
% Comparison of the performance of LMS and KLMS in nonlinear channel
% equalization. Execution time: < 1 minute (Intel Pentium Core2 Duo).
% 
% Weifeng Liu, Jose C. Principe and Simon Haykin, "Kernel  Adaptive
% Filtering: A Comprehensive Introduction", Wiley, 2010.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab
% https://github.com/steven2358/kafbox/

clear
close all

%% PARAMETERS
embedding = 5; % input embedding
delay = 2; % equalization delay

N_tr = 1000; % number of training data
N_te = 50; % number of testing data

setups{1} = lms(struct('mu',.01));
setups{2} = klms(struct('eta',0.2,'M',N_tr,'kerneltype','gauss','kernelpar',1));

%% PREPARE DATA

fprintf('Fig. 2.12 from "Kernel Adaptive Filtering: A Comprehensive\n');
fprintf('Introduction".\n')

u = randn(N_tr+N_te+embedding-1,1)>0;
u = 2*u-1; % binary input

z = u + 0.5*[0;u(1:end-1)]; % output of linear channel
ns = 0.4*randn(length(u),1); % channel noise
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
    fprintf('\n');
end

%% OUTPUT

figure
plot(MSE,'LineWidth',2)

legend('LMS','KLMS')
xlabel('iteration')
ylabel('testing MSE')
