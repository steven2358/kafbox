% Installation file. Adds local folders to path.

fprintf('Adding KAFBOX folders to Matlab path... ')

addpath(genpath([pwd filesep 'data']));
addpath(genpath([pwd filesep 'lib']));

addpath([pwd filesep 'demo']);

fprintf('done.\n')
disp('Type "savepath" if you wish to store changes.')
% savepath;
