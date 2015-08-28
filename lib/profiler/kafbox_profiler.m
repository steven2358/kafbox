% KAFBOX_PROFILER Profiler program. Runs through different algorithms, each
% with a number of configurations, each for a number of simulations.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function [data,algorithms,results] = ...
    kafbox_profiler(data,sim_opts,algorithms,outputdir)

num_algo = length(algorithms);
numsim = 1;
if isfield(sim_opts,'numsim'),
    numsim = sim_opts.numsim;
end

% pre-allocate structure
results = cell(num_algo,1);
for algo_ind = 1:num_algo,
    num_config = length(algorithms{algo_ind}.options.sweep_val);
    algo_results = cell(num_config,1);
    for config_ind = 1:num_config,
        algo_results{config_ind} = cell(numsim,1);
    end
    results{algo_ind} = algo_results;
end

for sim_ind = 1:numsim, % process simulations
    fprintf('SIM %2d\n',sim_ind);
    
    for algo_ind = 1:num_algo % process algorithms
        
        algo = algorithms{algo_ind};
        sweep_par = algo.options.sweep_par;
        sweep_val = algo.options.sweep_val;
        num_config = length(sweep_val);	% number of iterations in sweep
        
        for config_ind = 1:num_config, % process configurations
            algo_config = algo;
            fprintf('%9s %10s: ',algo_config.name,...
                struct2str(struct(sweep_par,sweep_val(config_ind))));
            
            % set sweep parameter value
            eval(sprintf('algo_config.options.%s = sweep_val(config_ind);',sweep_par));
            
            % perform one simulations for this configuration
            simresults = kafbox_profiler_simulation(data,sim_opts,...
                algo_config,sim_ind,outputdir);
            
            results{algo_ind}{config_ind}{sim_ind} = simresults;
            
            fprintf(' %.2fs',simresults.elapsed);
            fprintf('\n');
        end
    end
    fprintf('\n');
end
