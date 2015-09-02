% Installation file. Adds local folders to path.

fprintf('Adding KAFBOX folders to Matlab path... ')

addpath(genpath(fullfile(pwd,'data'))); % add data folder with subfolders
addpath(genpath(fullfile(pwd,'lib'))); % add lib folder with subfolders

addpath(fullfile(pwd,'demo')); % add demo folder without subfolders

fprintf('done.\n')
disp('Type "savepath" if you wish to store the changes.')
% savepath;
