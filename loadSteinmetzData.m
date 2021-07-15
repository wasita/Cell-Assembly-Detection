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

% Sessions given from colab ver of dataset
% CA1 sessions: 1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 18, 19, 22, 23, 26, 27, 29, 32, 34, 38
% CA3: 0, 6, 7, 8, 14, 15, 17, 30, 32, 37
% ACA sessions: 0, 3, 4, 11, 12, 21, 24, 26, 29, 34, 38

CA1_sessions = [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 18, 19, 22, 23, 26, 27, 29, 32, 34, 38];
CA3_sessions = [0, 6, 7, 8, 14, 15, 17, 30, 32, 37];
ACA_sessions = [0, 3, 4, 11, 12, 21, 24, 26, 29, 34, 38];

badSessions = [];

for session = 1:10
    
    try
        dataStruct = [];
        
        currFolder = steinmetz_data_dir(session+1);
        currSession = currFolder.name;
        currFolderContents = dir(fullfile([currFolder.folder, '/', currFolder.name], '*.npy'));
        
        fprintf(['Session #', num2str(session), ' ', currSession]);
        
        nNPYs = length(currFolderContents);
        
        % peakChannels is 1085 long of values from 1-740
        % each row is a neuron/cluster; value corresponds to cluster it's
        % coming from
        % should be same number of neurons!!! but it not :(
        peakChannelsFileIdx = strcmp({currFolderContents(:).name}, 'clusters.peakChannel.npy');
        peakChannels = readNPY([currFolderContents(peakChannelsFileIdx).name]); % 0; session 3 has 1085
        
        % number of channels varies by session -- session 3 has 1122
        % python says there should be 1769 for session 3
        channels = tdfread([steinmetz_data_path, '/', currSession , '/channels.brainLocation.tsv']);
        
        spikeClustersFileIdx =  strcmp({currFolderContents(:).name}, 'spikes.clusters.npy');
        spikeClusters = readNPY([currFolderContents(spikeClustersFileIdx).name]); % max val is 1084 or 1085
        
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
        
        maxSpikeTime = max(spikeTimes);
        binSize = round(maxSpikeTime/0.01); % bc times is given to us in seconds
        binRanges = 0:0.01:maxSpikeTime;
        
        
        %%% ACA
        
        if ismember(session, ACA_sessions)
            
            % number of channels of the probe corresponding to each area
            % session 3 channel indices range from 909 - 1058;
            % out of bounds of peakChannels
            ACA_channel_idx = find(strcmp(cellstr(channels.allen_ontology), 'ACA'));
            
            % number of neurons recorded from each area
            % session = 3 -- not finding ACA channel idx in peakchannels
            is_ACA = ismember(peakChannels, ACA_channel_idx);
            
            ACA_timesPerCluster = timesPerCluster(is_ACA);
            
            % bin into 10 ms bins
            ACA_mat = nan(length(ACA_timesPerCluster),binSize+1);
            
            for c = 1:length(ACA_timesPerCluster)
                currentCell = ACA_timesPerCluster(c);
                
                binCounts = histc([currentCell{:}],binRanges);
                ACA_mat(c,:) = binCounts;
            end
            
            [ACA_AssemblyTemplates, n_assemblies_ACA] = assembly_patterns(ACA_mat,opts);
            
            figure;
            subplot(211)
            stem(ACA_AssemblyTemplates(:,1))
            title([currSession, newline, 'ACA Assembly Patterns', newline, opts.threshold.method, ...
                ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
                num2str(n_assemblies_ACA)]);
            prettifyFig();
            subplot(212)
            stem(ACA_AssemblyTemplates(:,2))
            prettifyFig();
            
            saveas(gcf, [base_path, '/figures/', currSession, '_ACA_assembly_patterns_', ...
                opts.threshold.method, '_', opts.Patterns.method, '_', num2str(n_assemblies_ACA), ...
                '_assemblies.png']);
            
        end
        
        
        
        %%% CA3
        
        if ismember(session, CA3_sessions)
            
            % number of channels of the probe corresponding to each area
            CA3_channel_idx = find(strcmp(cellstr(channels.allen_ontology), 'CA3'));
            
            % number of neurons recorded from each area
            is_CA3 = ismember(peakChannels, CA3_channel_idx);
            
            CA3_timesPerCluster = timesPerCluster(is_CA3);
            
            % bin into 10 ms bins
            CA3_mat = nan(length(CA3_timesPerCluster),binSize+1);
            
            for c = 1:length(CA3_timesPerCluster)
                currentCell = CA3_timesPerCluster(c);
                
                binCounts = histc([currentCell{:}],binRanges);
                CA3_mat(c,:) = binCounts;
            end
            
            [CA3_AssemblyTemplates, n_assemblies_CA3] = assembly_patterns(CA3_mat,opts);
            
            % Plot detected assemblies
            figure;
            subplot(211)
            stem(CA3_AssemblyTemplates(:,1))
            title([currSession, newline, 'CA3 Assembly Patterns', newline, opts.threshold.method, ...
                ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
                num2str(n_assemblies_CA3)]);
            prettifyFig();
            subplot(212)
            stem(CA3_AssemblyTemplates(:,2))
            prettifyFig();
            
            saveas(gcf, [base_path, '/figures/', currSession, '_CA3_assembly_patterns_', ...
                opts.threshold.method, '_', opts.Patterns.method, '_', num2str(n_assemblies_CA3), ...
                '_assemblies.png']);
            
        end
        
        
        
        %%% CA1
        
        if ismember(session, CA1_sessions)
            
            % number of channels of the probe corresponding to each area
            CA1_channel_idx = find(strcmp(cellstr(channels.allen_ontology), 'CA1'));
            
            % number of neurons recorded from each area
            is_CA1 = ismember(peakChannels, CA1_channel_idx);
            
            CA1_timesPerCluster = timesPerCluster(is_CA1);
            
            % bin into 10 ms bins
            CA1_mat = nan(length(CA1_timesPerCluster),binSize+1);
            
            for c = 1:length(CA1_timesPerCluster)
                currentCell = CA1_timesPerCluster(c);
                
                binCounts = histc([currentCell{:}],binRanges);
                CA1_mat(c,:) = binCounts;
            end
            
            [CA1_AssemblyTemplates, n_assemblies_CA1] = assembly_patterns(CA1_mat,opts);
            
            % Plot detected assemblies
            figure;
            subplot(211)
            stem(CA1_AssemblyTemplates(:,1))
            title([currSession, newline, 'CA1 Assembly Patterns', newline, opts.threshold.method, ...
                ' with ', opts.Patterns.method, newline, 'Assemblies detected: ', ...
                num2str(n_assemblies_CA1)]);
            prettifyFig();
            subplot(212)
            stem(CA1_AssemblyTemplates(:,2))
            prettifyFig();
            
            saveas(gcf, [base_path, '/figures/', currSession, '_CA1_assembly_patterns_', ...
                opts.threshold.method, '_', opts.Patterns.method, '_', num2str(n_assemblies_CA1), ...
                '_assemblies.png']);
            
        end
        
    catch me
        
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
        
        badSessions = [badSessions, session];
        
    end
    
end


