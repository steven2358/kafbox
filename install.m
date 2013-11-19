% Installation file. Adds local folders to path.

fprintf('Adding KAFBOX folders to Matlab path... ')

addpath(genpath(fullfile(pwd,'data')));
addpath(genpath(fullfile(pwd,'lib')));

addpath((fullfile(pwd,'demo')));

fprintf('done.\n')
disp('Type "savepath" if you wish to store changes.')
% savepath;
