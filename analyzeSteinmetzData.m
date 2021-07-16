%% loadSteinmetzData

%%% TO DO %%%
% [] Try different time bins; iterate and step through different bins
% [] Try ICA over PCA
% [] Make things task-related (Looking at correct/incorrect trials; task
% segments)
% [] PSCTHs -- correct correlate to choice-time


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


%% Set up which sessions to iterate over
% Sessions given from colab ver of dataset
CA1_sessions = [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 18, 19, 22, 23, 26, 27, 29, 32, 34, 38];
CA3_sessions = [0, 6, 7, 8, 14, 15, 17, 30, 32, 37];
ACA_sessions = [0, 3, 4, 11, 12, 21, 24, 26, 29, 34, 38];

% force 0-index to match python colab data indices
allSessions = 0:nFolders-1;

% check if we have a badSessions already
% badSessions = [3     4     7    11    12    13    14    15    29    32    34    38]
% if yes, load it
% if no, maybe assume it needs to be re-run?
% and also do not iterate over those sessions in for loop
% if exist([base_path, '/badSessions.mat'], 'file') ~= 0
%     loadedBadSessions = load([base_path, '/badSessions.mat']);
%     badSessions = loadedBadSessions.badSessions;
%     goodSessions = allSessions(~ismember(allSessions,badSessions)); % remove badSessions from sessions to iterate over
% else
%     badSessions = [];
% end

%% Pattern assembly extraction parameters %%

methods = {'PCA', 'ICA'};
timeWindows = [0.25:.25:1; ];

opts.threshold.permutations_percentile = 95;
opts.threshold.number_of_permutations = 20;
opts.threshold.method = 'MarcenkoPastur';

%% Loop through each subfolder in steinmetz-data and load in .npy in each folder
for session = 0:nFolders-1
    
    % Iterate over different methods
    for method = 1:length(methods)
        
        patternMethod = methods{method};
        
        try
            dataStruct = [];
            
            currFolderDir = steinmetz_data_dir(session+1);
            currSession = currFolderDir.name;
            sessionPath = [currFolderDir.folder, '/', currFolderDir.name];
            
            fprintf(['Session #', num2str(session), ' ', currSession, newline]);
            
            nNPYs = length(currFolderContents);
            
            % read in cluster annotation file
            %  1 = MUA -- presumed to contain spikes from multiple analyzed,
            %  not analyzed
            % 2 = good
            % 3 = unsorted
            clusterAnnotes = readNPY([sessionPath, '/clusters._phy_annotation.npy']);
            
            % each row is a neuron/cluster; value corresponds to cluster it's
            % coming from
            peakChannels = readNPY([sessionPath, '/clusters.peakChannel.npy']);
            
            % Sahana's way to confirm we get same result for peakChannels:
            %             peakChannels = readNPY(sprintf( '%s/%s/clusters.peakChannel.npy',steinmetz_data_path, currFolder.name));
            
            % number of channels varies by session
            channels = tdfread([sessionPath , '/channels.brainLocation.tsv']);
            spikeClusters = readNPY([sessionPath, '/spikes.clusters.npy']);
            spikeTimes =readNPY([sessionPath, '/spikes.times.npy']);
            
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
            binStep = 0.1; % first tried 0.01 for 10ms; now we try 0.1 = 100 ms
            binRanges = 0:binStep:maxSpikeTime;
            
            % Set parameters for running assemply_patterns
            % Moved these 3 to outside of loop as they're constants
            %             opts.threshold.permutations_percentile = 95;
            %             opts.threshold.number_of_permutations = 20;
            %             opts.threshold.method = 'MarcenkoPastur';
            opts.Patterns.method = patternMethod;
            
            
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
                    num2str(n_assemblies_ACA)], newline, 'Bin Step = ', num2str(binStep));
                prettifyFig();
                subplot(212)
                stem(ACA_AssemblyTemplates(:,2))
                prettifyFig();
                set(gcf, 'Visible', 'off');
                
                saveas(gcf, [base_path, '/figures/AssemblyPatterns/', currSession, '_ACA_assembly_patterns_', ...
                    opts.threshold.method, '_', opts.Patterns.method, '_', 'binStep_', ...
                    num2str(binStep), '_', num2str(n_assemblies_ACA), ...
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
                    num2str(n_assemblies_CA3), newline, 'Bin Step = ', num2str(binStep)]);
                prettifyFig();
                subplot(212)
                stem(CA3_AssemblyTemplates(:,2))
                prettifyFig();
                set(gcf, 'Visible', 'off');
                
                saveas(gcf, [base_path, '/figures/AssemblyPatterns/', currSession, '_CA3_assembly_patterns_', ...
                    opts.threshold.method, '_', opts.Patterns.method, '_', 'binStep_', ...
                    num2str(binStep), '_', num2str(n_assemblies_CA3), ...
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
                    num2str(n_assemblies_CA1), newline, 'Bin Step = ', num2str(binStep)]);
                prettifyFig();
                subplot(212)
                stem(CA1_AssemblyTemplates(:,2))
                prettifyFig();
                set(gcf, 'Visible', 'off');
                
                saveas(gcf, [base_path, '/figures/AssemblyPatterns/', currSession, '_CA1_assembly_patterns_', ...
                    opts.threshold.method, '_', opts.Patterns.method, '_', 'binStep_', ...
                    num2str(binStep), '_', num2str(n_assemblies_CA1), ...
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
    
end

%% Save out which sessions are bad to load in later if needed
% to not need to re-run the above for loop
save('badSessions.mat', 'badSessions');
