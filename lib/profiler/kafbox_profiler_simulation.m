% Perform a single profiler simulation.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

function simresults = kafbox_profiler_simulation(data,sim_opts,...
    algo_config,sim_ind,output_dir)

simdata = data;

if isfield(data,'offset')
    simdata.offset = (sim_ind-1)*data.offset;
elseif isfield(data,'permutation')
    simdata.permutation = (sim_ind-1)*data.permutation; % start with 0 = no permutation
else
    % fix random seed per simulation
    rs = sim_ind;
    rng('default');
    rng(sim_ind);
    simdata.rs = rs;
end

% check if results for this simulation have been stored before
simresults = kafbox_profiler_storet(simdata,algo_config,output_dir);

if isempty(simresults), % perform simulation

    [X,Y,X_test,Y_test] = kafbox_data(simdata); % load data
    
    every = 1;
    if isfield(data,'test_every')
        every = data.test_every;
    end
    
    N = size(X,1);
    all_fl = zeros(N,1); %  flops per iteration
    all_bytes = zeros(N,1); %  flops per iteration
    erm = sim_opts.error_measure; % error measure
    eval(sprintf('%s = nan*zeros(N,1);',erm));
    var_test = var(Y_test); % for NMSE
    
    kaf = feval(sprintf('%s_profiler',algo_config.class),algo_config.options);
    
    results2store = {'elapsed','flops','bytes',erm};

    y_test_ind = ones(N,1);
    if isfield(data,'N_switch')
        y_test_ind(data.N_switch:end) = 2;
    end
    
    for i=1:N,
        Y_est = kaf.evaluate(X(i,:)); % predict
        kaf = kaf.train_profiled(X(i,:),Y(i)); % train
        
        all_fl(i) = kaf.lastflops();
        all_bytes(i) = kaf.lastbytes();
        
        if isfield(data,'test_every_conv')
            de = data.test_every_conv;
            if i<10*data.test_every_conv
                % start from 1 at i=1 and go exponentially to every at
                % i=10*every
                every = max(1,round(10^(log10(de)/10/de*i)));
            else
                every = data.test_every_conv;
            end
        end
        if mod(i,every) == 0
            if isempty(X_test)
                Y_test = Y(i); % test prediction
            else
                Y_est = kaf.evaluate(X_test); % test on multiple data
            end
            erm_val = mean((Y_test(:,y_test_ind(i)) - Y_est).^2);
            if strcmp(erm,'NMSE'),
                erm_val = erm_val/var_test; %#ok<NASGU>
            end
            eval(sprintf('%s(i) = erm_val;',erm))
        end
    end
    
    elapsed = kaf.elapsed; %#ok<NASGU>
    flops = max(all_fl); %#ok<NASGU>
    bytes = max(all_bytes); %#ok<NASGU>
    
    for j=1:length(results2store)
        eval(sprintf('simresults.%s = %s;',results2store{j},results2store{j}))
    end
    
    kafbox_profiler_storet(simdata,algo_config,output_dir,simresults);
    
    fprintf('calculated');
else
    fprintf('retrieved ');
end
