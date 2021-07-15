%% loadSteinmetzData

%% housekeeping
clear all; clc;


%% Set up directories 

%%% CHANGE YOUR BASE_PATH HERE %%% 
% set up main dir
base_path = '/Users/wasita/Documents/GitHub/Cell-Assembly-Detection';
steinmetz_data_path = [base_path, '/steinmetz-data'];

addpath(genpath(base_path))

% add npy-matlab to path
% found that genpath worked better than the method suggested by their
% readme
addpath(genpath([base_path,'/npy-matlab']));

% set up steinmetz-data dir
steinmetz_data_dir = dir(steinmetz_data_path);

% if you're on a PC, can change to indexing up to 2 instead of 3
steinmetz_data_dir(1:3) = []; % ignore current, parent dirs and .DS_Store
nFolders = length(steinmetz_data_dir);

%% Loop through each subfolder in steinmetz-data and load in .npy in each folder

for session = 1:nFolders
    
    dataStruct = [];
    
    currFolder = steinmetz_data_dir(session);
    currSession = currFolder.name;
    currFolderContents = dir(fullfile([currFolder.folder, '/', currFolder.name], '*.npy'));
    
    nNPYs = length(currFolderContents);

    peakChannelsFileIdx = strcmp({currFolderContents(:).name}, 'clusters.peakChannel.npy');
    peakChannels = readNPY([currFolderContents(peakChannelsFileIdx).name]);
    
    % NOTE: currently channels is specific to 1 session but the outer loop
    % is iterating over all sessions
    channels = tdfread([steinmetz_data_path, '/', currSession , '/channels.brainLocation.tsv']);
   
    spikeClustersFileIdx =  strcmp({currFolderContents(:).name}, 'spikes.clusters.npy');
    spikeClusters = readNPY([currFolderContents(spikeClustersFileIdx).name]);
    
    spikeTimesFileIdx = strcmp({currFolderContents(:).name}, 'spikes.times.npy');
    spikeTimes = readNPY([currFolderContents(spikeTimesFileIdx).name]);
    
    % add spikeClusters and spikeTimes to dataStruct to later turn into a
    % data table
    dataStruct.spikeClusters = spikeClusters;
    dataStruct.spikeTimes = spikeTimes;
    
    dataTable = struct2table(dataStruct);
    
    % get number of peak channels
    nPeakChannels = max(peakChannels);
    
    % sort table
    sortedTable = sortrows(dataTable, 'spikeClusters');
    
    % grab all the spike times for a given cluster
    % the size of this varies so we'll initialize without pre-allocation
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
    % and plotting results of found assemblies
    
    % set up parameters for MarcenkoPastur PCA for ACA & CA3
    opts.threshold.permutations_percentile = 95;
    opts.threshold.number_of_permutations = 20;
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'PCA';
    
    [ACA_AssemblyTemplates, n_assemblies_ACA] = assembly_patterns(ACA_mat,opts);
    [CA3_AssemblyTemplates, n_assemblies_CA3] = assembly_patterns(CA3_mat,opts);
    
    figure;
    subplot(211)
    stem(ACA_AssemblyTemplates(:,1))
    title(['ACA Assembly Patterns', newline, opts.threshold.method, ...
        ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
        num2str(n_assemblies_ACA)]);
    prettifyFig();
    subplot(212)
    stem(ACA_AssemblyTemplates(:,2))
    prettifyFig();
    
    figure;
    subplot(211)
    stem(CA3_AssemblyTemplates(:,1))
    title(['CA3 Assembly Patterns', newline, opts.threshold.method, ...
        ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
        num2str(n_assemblies_CA3)]);
    prettifyFig();
    subplot(212)
    stem(CA3_AssemblyTemplates(:,2))
    prettifyFig();
    
    
    % set up parameters for circularshift PCA for ACA & CA3
    opts.threshold.permutations_percentile = 95;
    opts.threshold.number_of_permutations = 20;
    opts.threshold.method = 'circularshift';
    opts.Patterns.method = 'PCA';
    
    [ACA_AssemblyTemplates, n_assemblies_ACA] = assembly_patterns(ACA_mat,opts);
    [CA3_AssemblyTemplates, n_assemblies_CA3] = assembly_patterns(CA3_mat,opts);
    
    figure;
    subplot(211)
    stem(ACA_AssemblyTemplates(:,1))
    title(['ACA Assembly Patterns', newline, opts.threshold.method, ...
        ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
        num2str(n_assemblies_ACA)]);
    prettifyFig();
    subplot(212)
    stem(ACA_AssemblyTemplates(:,2))
    prettifyFig();
    
    figure;
    subplot(211)
    stem(CA3_AssemblyTemplates(:,1))
    title(['CA3 Assembly Patterns', newline, opts.threshold.method, ...
        ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
        num2str(n_assemblies_CA3)]);
    prettifyFig();
    subplot(212)
    stem(CA3_AssemblyTemplates(:,2))
    prettifyFig();
    
    
    % TO ADD BACK IN LATER
    %     Activities = assembly_activity(AssemblyTemplates,ACA_mat);
    %
    %     figure(4),clf
    %     subplot(211)
    %     imagesc(ACA_mat)
    %     xlim([0 100])
    %     subplot(212)
    %     plot(Activities')
    %     xlim([0 100])

    
end


