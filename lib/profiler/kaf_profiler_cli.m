% Profiler program.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function [results,setups] = kaf_profiler_cli(config_file,output_dir) %#ok<STOUT>
% KAF_PROFILER_CLI Summary
%   Details

fprintf(sprintf('Sweep on %s configuration.\n',config_file));
eval(sprintf('%s',config_file));

if nargin < 2
    savedata = 0;
else
    savedata = 1;
end


num_setup = length(setups);
% titles = cell(num_setup,1);
results.flops = cell(num_setup,1);
results.bytes = cell(num_setup,1);
results.msedb = cell(num_setup,1);

for setup_ind = 1:length(setups)
    
    setup = setups{setup_ind};
    %     titles{setup_ind} = setup.algo.name;
    sw_par = setup.algo.options.sw_par;
    sw_val = setup.algo.options.sw_val;
    num_sw = length(sw_val);	% number of iterations in sweep
    
    flops = zeros(num_sw,1);
    bytes = zeros(num_sw,1);
    msedb = zeros(num_sw,1);
    
    for sw_ind = 1:num_sw,
        fprintf('%s %s=%.3f ',setup.algo.name,sw_par,sw_val(sw_ind));
        % pd = setup.data;
        
        eval(sprintf('setup.algo.options.%s = sw_val(sw_ind);',sw_par));
        
        % check if results for this setup have been stored before
        setupresults = kaf_setuphandler(setup,output_dir);
        
        if isempty(setupresults), % perform simulation
            [X,Y] = kaf_getdata(setup.data); % load data
            
            N = size(X,1);
            
            all_fl = zeros(N,1); %  flops per iteration
            all_bytes = zeros(N,1); %  flops per iteration
            SE = zeros(N,1);
            
            kaf = feval(sprintf('%s_profiler',setup.algo.name),setup.algo.options);
            
            t1 = tic;
            for i=1:N,
                if ~mod(i,floor(N/10)), fprintf('.'); end
                y_est = kaf.evaluate(X(i,:)); % evaluate
                kaf = kaf.train(X(i,:),Y(i)); % train
                
                all_fl(i) = kaf.lastflops();
                all_bytes(i) = kaf.lastbytes();
                
                SE(i) = (Y(i)-y_est)^2; % calcalate test error
            end
            
            setupresults.flops = mean(all_fl);
            setupresults.bytes = max(all_bytes);
            setupresults.msedb = 10*log10(mean(SE(1000:end)));
            
            kaf_setuphandler(setup,output_dir,setupresults);
            
            t2 = toc(t1);
            fprintf(1,' %.2fs.' ,t2);
        else
            fprintf(' _______');
        end
        
        flops(sw_ind) = setupresults.flops;
        bytes(sw_ind) = setupresults.bytes;
        msedb(sw_ind) = setupresults.msedb;
        
        fprintf(1,' Average:');
        fprintf(1,' %.2f KFLOPS',flops(sw_ind)/1000);
        fprintf(1,' %d bytes max',bytes(sw_ind));
        
        fprintf(1,', MSE=%.2fdB\n',msedb(sw_ind));
        
    end
    fprintf('\n');
    results.flops{setup_ind} = flops;
    results.bytes{setup_ind} = bytes;
    results.msedb{setup_ind} = msedb;
end

