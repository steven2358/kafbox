% Handles storage and retrieval of profiler results.
%
% Checks if the results for a given configuration were stored before, and
% saves new results if necessary.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/ 

function results = kafbox_profiler_storet(data,config,output_dir,results)

if nargin<4
    option = 'check';
else
    option = 'save';
end

dataset = lower(data.name);
if isfield(data,'class')
    dataset = sprintf('%s_%s',lower(data.class),dataset);
end
dataset_algo = [dataset '_' config.class];
results_path = fullfile(output_dir,dataset_algo);
index_path = fullfile(results_path,'index.mat');
if ~exist(results_path,'file')
    mkdir(results_path);
end

% all relevant values in one array
my_config = config.options;
fn = fieldnames(data);
for i = 1:length(fn)
    if ~strcmp(fn{i},'name')
        my_config.(sprintf('data_%s',fn{i})) = data.(fn{i});
    end
end

skiplist = {'sweep_par','sweep_val','data_numsim'};
for i = 1:length(skiplist)
    if isfield(my_config,skiplist{i})
        my_config = rmfield(my_config,skiplist{i});
    end
end

switch option
    case 'check'
        if exist(index_path,'file')==2
            load(index_path); % loads "configs"
            results = find_results(my_config,configs,results_path); %#ok<NODEF>
        else
            results = [];
        end
    case 'save'
        if exist(index_path,'file')==2
            load(index_path); % loads "configs"
            ind = length(configs)+1; %#ok<NODEF>
        else
            ind = 1;
        end
        id = datestr(now,30);
        configs{ind}.id = id;
        configs{ind}.cstr = struct2str(my_config);
        save(index_path,'configs');
        save_results(results,results_path,id)
    otherwise
        error('unknown option')
end


% check if config was already processed and return corresponding results
function r = find_results(my_config,configs,results_path)

str = struct2str(my_config);

for i=1:length(configs)
    stri = configs{i}.cstr;
    if strcmp(str,stri)
        id = configs{i}.id;
        r = load_results(results_path,id);
        return
    end
end
r = [];


function save_results(results,results_path,id)
fname = sprintf('%s/%s',results_path,id);
save(fname,'results');


function results = load_results(results_path,id)
fname = sprintf('%s/%s.mat',results_path,id);
if exist(fname,'file')==2
    r = load(fname);
    results = r.results;
else
    results = [];
end
