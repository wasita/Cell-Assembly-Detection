%% loadSteinmetzData

%% Set up directories 

% set up main dir
addpath(genpath('/Users/wasita/Documents/GitHub/Cell-Assembly-Detection'))

% add npy-matlab to path
addpath(genpath('/Users/wasita/Documents/GitHub/Cell-Assembly-Detection/npy-matlab'))

% set up steinmetz-data dir
steinmetz_data_dir = dir('steinmetz-data');
steinmetz_data_dir(1:3) = []; % ignore current, parent dirs and .DS_Store % if you're on a PC, can change to indexing up to 2 instead of 3
nFolders = length(steinmetz_data_dir);

%% Loop through each subfolder and load in .npy in each folder

for folder = 1:nFolders
    
    currFolder = steinmetz_data_dir(folder);
    currFolderContents = dir(fullfile([currFolder.folder, '/', currFolder.name], '*.npy'));
    
    nNPYs = length(currFolderContents);
    
    for n = 1:nNPYs
        
        loadedNPY = readNPY([currFolderContents(n).name]);
        
    end
    
    
end