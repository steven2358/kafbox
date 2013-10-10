% Checks if the results for a given configuration were stored before, and
% saves new results if necessary.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function results = kafbox_resultshandler6(data,config,output_dir,results)
% store all results for one algorithm in a single file

if nargin<4
    option = 'check';
else
    option = 'save';
end

if isfield(data,'file')
    [pathstr,dataset] = fileparts(data.file); %#ok<ASGLU>
elseif isfield(setup.data,'generate')
    dataset = data.generate;
end
index_path = sprintf('%s/%s_%s_index.mat',output_dir,dataset,config.class);
results_path = sprintf('%s/results/%s_%s',output_dir,dataset,config.class);

if exist('output','file')
    if ~exist('output/results','file')
        mkdir('output/results');
        mkdir(results_path);
    else
        if ~exist(results_path,'file')
            mkdir(results_path);
        end
    end
else
    mkdir('output');
    mkdir('output/results');
    mkdir(results_path);
end


% all relevant values in one array
my_config = config.options;
fn = fieldnames(data);
for i = 1:length(fn),
    if ~strcmp(fn{i},'name')
        my_config.(sprintf('data_%s',fn{i})) = data.(fn{i});
    end
end

skiplist = {'sweep_par','sweep_val','data_numsim'};
for i = 1:length(skiplist),
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
        configs{ind}.cstr = struct2str(my_config); %#ok<NASGU>
        save(index_path,'configs');
        save_results(results,results_path,id)
    otherwise
        error('unknown option')
end



% check if config was already processed and return corresponding results
function r = find_results(my_config,configs,results_path)

str = struct2str(my_config);

for i=1:length(configs),
    stri = configs{i}.cstr;
    if strcmp(str,stri),
        id = configs{i}.id;
        r = load_results(results_path,id);
        return
    end
end
r = [];



function save_results(results,results_path,id) %#ok<INUSL>
fname = sprintf('%s/%s',results_path,id);
save(fname,'results');



function results = load_results(results_path,id)
fname = sprintf('%s/%s',results_path,id);
r = load(fname);
results = r.results;



% convert a structure to a string
function str = struct2str(my_struct)
fields = sortrows(fieldnames(my_struct)); % avoid permutations
str = '';
for i=1:length(fields)
    fn = fields{i};
    fv = my_struct.(fn);
    switch class(fv)
        case 'char'
            str = sprintf('%s_%s=%s',str,fn,fv);
        case 'double'
            if (round(fv)==fv)
                str = sprintf('%s_%s=%d',str,fn,fv);
            else
                found = 0;
                for j=0:10,
                    if (round(10^j*fv)==10^j*fv) && ~found
                        str = sprintf(sprintf('%%s %%s=%%.%df',j),str,fn,fv);
                        found = 1;
                    end
                end
            end
        otherwise
            error('unknown field class');
    end
end
str = str(2:end);

