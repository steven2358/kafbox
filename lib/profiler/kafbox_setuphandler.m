% Checks if the results for a given setup were stored before, and saves new
% results if necessary.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function setupresults = kafbox_setuphandler(setup,output_dir,setupresults)
% store all results for one algorithm in a single file

if nargin<3
    option = 'check';
else
    option = 'save';
end

if isfield(setup.data,'file')
    [pathstr,dataset] = fileparts(setup.data.file); %#ok<ASGLU>
elseif isfield(setup.data,'generate')
    dataset = setup.data.generate;
end
fname = sprintf('%s/%s_%s.mat',output_dir,dataset,setup.algo.name);

% all relevant values in one array
my_setup = setup.algo.options;
fn = fieldnames(setup.data);
for i = 1:length(fn),
    if ~strcmp(fn{i},'name')
        my_setup.(sprintf('data_%s',fn{i})) = setup.data.(fn{i});
    end
end

skiplist = {'sweep_par','sweep_val'};
for i=1:length(skiplist),
    my_setup = rmfield(my_setup,skiplist{i});
end

switch option
    case 'check'
        if exist(fname,'file')==2
            load(fname); % loads the variables setups and results
            setupresults = find_setupresults(my_setup,setups,results); %#ok<NODEF>
        else
            setupresults = [];
        end
    case 'save'
        if exist(fname,'file')==2
            load(fname); % loads the variables setups and results
            ind = length(setups)+1; %#ok<NODEF>
        else
            ind = 1;
        end
        setups{ind} = my_setup; %#ok<NASGU>
        results{ind} = setupresults; %#ok<NASGU>
        save(fname,'setups','results');
    otherwise
        error('unknown option')
end


% check if setup was already processed and return corresponding results
function r = find_setupresults(setup,setups,results)

for i=1:length(setups),
    setupi = setups{i};
    
    fieldsi = sortrows(fieldnames(setupi));
    fields = sortrows(fieldnames(setup));
    identical = true;

    for f = 1:length(fields)
        if isfield(setupi,fields{f})
            field = sprintf('setup.%s',fields{f});
            fieldi = sprintf('setupi.%s',fields{f});
            if eval(field) ~= eval(fieldi)
                identical = false;
            end
        else
            identical = false;
        end
    end
    
    for fi = 1:length(fieldsi)
        if isfield(setup,fieldsi{fi})
            field = sprintf('setup.%s',fieldsi{fi});
            fieldi = sprintf('setupi.%s',fieldsi{fi});
            if eval(field) ~= eval(fieldi)
                identical = false;
            end
        else
            identical = false;
        end
    end
    
    if identical,
        r = results{i};
        return
    end
end
r = [];
