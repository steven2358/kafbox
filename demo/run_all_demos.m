% Script to run all demos consecutively.
%
% Author: Steven Van Vaerenbergh, 2015.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox

clc
close all
clear

% get list of demo scripts
fdir = fileparts(which('run_all_demos.m'));
files_demo = dir(fullfile(fdir,'demo_*.m'));
files = files_demo;
folders = repmat({fdir},length(files),1);

% get list of literature scripts
folders_literature = dir(fullfile(fdir,'literature','*2*'));
for i=1:length(folders_literature)
    folder_i = fullfile(fdir,'literature',folders_literature(i).name);
    files_i = dir(fullfile(folder_i,'fig*.m'));
    folders_i = repmat({folder_i},length(files_i),1);
    
    files = [files; files_i]; %#ok<AGROW>
    folders = [folders; folders_i]; %#ok<AGROW>
end

[~,files] = cellfun(@fileparts, {files.name}, 'UniformOutput',0);

t1 = tic;
fprintf('\n')
for i=1:length(files)
    close all
    clear eval
    save(fullfile(tempdir,'temp.mat'),'i','folders','files','fdir',...
        't1'); % memory map
    
    try
        % run script
        cd(folders{i})
        fname_demo = files{i};
        fprintf('\nRunning %s\n',fname_demo);
        eval(fname_demo);
    catch err
        % return to demo folder
        load(fullfile(tempdir,'temp.mat'));
        cd(fdir);
        
        % report error
        me = err.stack(1);
        error(['Error in ',...
            '<a href="matlab: opentoline(''%s'',%d)">',...
            '%s at %d</a>\n',...
            '%s'],...
           me.file, me.line, me.name, me.line,...
            err.message);
    end

    load(fullfile(tempdir,'temp.mat'));
    cd(fdir);
end
delete(fullfile(tempdir,'temp.mat'));
toc(t1)

close all

fprintf('\nAll tests finished.\n')
