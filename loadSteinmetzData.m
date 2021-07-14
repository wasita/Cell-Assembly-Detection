%% loadSteinmetzData

%% Set up directories 

% set up main dir
addpath(genpath('/Users/wasita/Documents/GitHub/Cell-Assembly-Detection'))

% add npy-matlab to path
% found that genpath worked better than the method suggested by their
% readme
addpath(genpath('/Users/wasita/Documents/GitHub/Cell-Assembly-Detection/npy-matlab'))

% set up steinmetz-data dir
steinmetz_data_dir = dir('steinmetz-data');
steinmetz_data_dir(1:3) = []; % ignore current, parent dirs and .DS_Store % if you're on a PC, can change to indexing up to 2 instead of 3
nFolders = length(steinmetz_data_dir);

%% Loop through each subfolder and load in .npy in each folder

dataStruct = [];

for folder = 1:nFolders
    
    currFolder = steinmetz_data_dir(folder);
    currFolderContents = dir(fullfile([currFolder.folder, '/', currFolder.name], '*.npy'));
    
    nNPYs = length(currFolderContents);

    peakChannels = readNPY([currFolderContents(10).name]);
    
    % NOTE: currently channels is specific to 1 session but the outer loop
    % is iterating over all sessions
    channels = tdfread('/Users/wasita/Documents/GitHub/Cell-Assembly-Detection/steinmetz-data/Cori_2016-12-14/channels.brainLocation.tsv');
   
    spikeClusters = readNPY([currFolderContents(33).name]);
    spikeTimes = readNPY([currFolderContents(35).name]);
    
    dataStruct.spikeClusters = spikeClusters;
    dataStruct.spikeTimes = spikeTimes;
    
    dataTable = struct2table(dataStruct);
    
    nPeakChannels = max(peakChannels);
    
    sortedTable = sortrows(dataTable, 'spikeClusters');
    
    % grab all the spike times for a given cluster
    %     timesPerCluster = cell(max(spikeClusters),1);
    timesPerCluster = {};
        
    for cluster_idx = 0:max(spikeClusters)
        
        allTimesForCurrCluster = spikeTimes(spikeClusters==cluster_idx);
        timesPerCluster{cluster_idx+1} = allTimesForCurrCluster;
       
    end

    % number of channels of the probe corresponding to each area
    ACA_channel_idx = find(strcmp(cellstr(channels.allen_ontology), 'ACA'));
    CA3_channel_idx = find(strcmp(cellstr(channels.allen_ontology), 'CA3'));
    
    % number of neurons recorded from each area
    is_ACA = ismember(peakChannels, ACA_channel_idx);
    is_CA3 = ismember(peakChannels, CA3_channel_idx);
    
    ACA_timesPerCluster = timesPerCluster(is_ACA);
    CA3_timesPerCluster = timesPerCluster(is_CA3);
    
    maxSpikeTime = max(spikeTimes);
    binSize = round(maxSpikeTime/0.01); % bc times is given to us in seconds
    binRanges = 0:0.01:maxSpikeTime;
    
    % bin ACA_timesPerCluster and CA3_timesPerCluster into 10 ms bins
    
    ACA_mat = nan(length(ACA_timesPerCluster),binSize+1);
    CA3_mat = nan(length(CA3_timesPerCluster),binSize+1);
    
    for c = 1:length(ACA_timesPerCluster)
        currentCell = ACA_timesPerCluster(c);
        
        binCounts = histc([currentCell{:}],binRanges);
        ACA_mat(c,:) = binCounts;
        
    end
    
    for c = 1:length(CA3_timesPerCluster)
        currentCell = CA3_timesPerCluster(c);
        
        binCounts = histc([currentCell{:}],binRanges);
        CA3_mat(c,:) = binCounts;
        
    end
    
    
    % try passing our matrices through assembly_patterns
    
    opts.threshold.permutations_percentile = 95;
    opts.threshold.number_of_permutations = 20;
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'PCA';
    AssemblyTemplates = assembly_patterns(ACA_mat,opts);
    
    
    opts.threshold.permutations_percentile = 95;
    opts.threshold.number_of_permutations = 20;
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'PCA';
    AssemblyTemplates = assembly_patterns(CA3_mat,opts);
    
    figure(2),clf
    subplot(211)
    stem(AssemblyTemplates(:,1))
    subplot(212)
    stem(AssemblyTemplates(:,2))

   
    Activities = assembly_activity(AssemblyTemplates,ACA_mat);
    
    figure(4),clf
    subplot(211)
    imagesc(ACA_mat)
    xlim([0 100])
    subplot(212)
    plot(Activities')
    xlim([0 100])

    
end


