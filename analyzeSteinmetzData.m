%% analyzeSteinmetzData

%%% TO DO %%%
% [x] Try different time bins; iterate and step through different bins

% [x] Try ICA over PCA

% [] Make things task-related (Looking at correct/incorrect trials; task
% segments)

% [] PSCTHs -- correct correlate to choice-time

% [x] Plot all assemblies not just first 2!

% [] Generate table for easier readout. Columns should be:
%   [] # assemblies detected
%   [] area (ACA, CA1, CA3)
%   [] method (PCA, ICA)
%   [] threshold method (MP)
%   [] time bins

% [] Activity plots

% [] Save out timesPerCluster per area per session

% [X] Refactor code -- make function to detect clusters per area

% [] Auto-generate directories for figures output

%% Housekeeping
clear all; clc;


%% Set up directories

%%% CHANGE YOUR BASE_PATH HERE %%%
% Set up main dir
base_path = '/Users/wasita/Documents/GitHub/Cell-Assembly-Detection';
steinmetz_data_path = [base_path, '/steinmetz-data'];

addpath(genpath(base_path))

% Add npy-matlab to path
% found that genpath worked better than the method suggested by their
% readme
addpath(genpath([base_path,'/npy-matlab']));

% Set up steinmetz-data dir
steinmetz_data_dir = dir(steinmetz_data_path);

%% NOTE: if you're on a PC, can change to indexing up to 2 instead of 3
% .DS_Store is a mac specific hidden file
steinmetz_data_dir(1:3) = []; % ignore current, parent dirs and .DS_Store
nFolders = length(steinmetz_data_dir);

%% Set up which sessions to iterate over
% Sessions given from colab ver of dataset
CA1_sessions = [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 18, 19, 22, 23, 26, 27, 29, 32, 34, 38];
CA3_sessions = [0, 6, 7, 8, 14, 15, 17, 30, 32, 37];
ACA_sessions = [0, 3, 4, 11, 12, 21, 24, 26, 29, 34, 38];

% Force 0-index to match python colab data indices
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
nMeths = length(methods);

ROIs = {'ACA', 'CA1', 'CA3'}; % should match label format in channels.ontology
nROIs = length(ROIs);

timeWindows = [0.025:0.025:.1, 0.25, 0.5];

opts.threshold.permutations_percentile = 95;
opts.threshold.number_of_permutations = 20;
opts.threshold.method = 'MarcenkoPastur';

%% Initialize table to read out results from batch running
outputTable = [];

%% Loop through each subfolder in steinmetz-data and load in .npy in each folder
for sessionIdx = 0:nFolders-1
    
    % Iterate over different methods
    for method = 1:nMeths
        
        patternMethod = methods{method};
        
        % We first tried 0.01 for 10ms; now we try 0.1 = 100 ms and stuff
        % in between
        for binStep = timeWindows
            
            try
                dataStruct = [];
                
                folderDir = steinmetz_data_dir(sessionIdx+1);
                sessionName = folderDir.name;
                
                % Example sessionPath:
                % '/Users/wasita/Documents/GitHub/Cell-Assembly-Detection/steinmetz-data/Cori_2016-12-14'
                sessionPath = [folderDir.folder, '/', folderDir.name];
                
                fprintf(['Session #', num2str(sessionIdx), ' ', sessionName, newline]);
                
                % Read in cluster annotation file
                %  1 = MUA -- presumed to contain spikes from multiple analyzed,
                %  not analyzed
                % 2 = good
                % 3 = unsorted
                clusterAnnotes = readNPY([sessionPath, '/clusters._phy_annotation.npy']);
                
                % Each row is a neuron/cluster; value corresponds to cluster it's
                % coming from
                % Number of peak channels (1st dim of peakChannels) varies
                % across sessions
                peakChannels = readNPY([sessionPath, '/clusters.peakChannel.npy']);

                % Number of channels varies by session
                channels = tdfread([sessionPath , '/channels.brainLocation.tsv']);
                spikeClusters = readNPY([sessionPath, '/spikes.clusters.npy']);
                spikeTimes =readNPY([sessionPath, '/spikes.times.npy']);
                
                % Add spikeClusters and spikeTimes to dataStruct to later turn into a
                % data table
                dataStruct.spikeClusters = spikeClusters;
                dataStruct.spikeTimes = spikeTimes;
                
                dataTable = struct2table(dataStruct);
                
                % Get number of peak channels
                nPeakChannels = max(peakChannels);
                
                % Sort table
                sortedTable = sortrows(dataTable, 'spikeClusters');
                
                % Grab all the spike times for a given cluster
                % the size of this varies so we'll initialize without pre-allocation
                timesPerCluster = {};
                
                % Number of clusters corresponds to length of 1st dim in
                % peakChannels
                for cluster_idx = 0:max(spikeClusters)
                    allTimesForCurrCluster = spikeTimes(spikeClusters==cluster_idx);
                    timesPerCluster{cluster_idx+1} = allTimesForCurrCluster;
                end
                
                maxSpikeTime = max(spikeTimes);
                binSize = round(maxSpikeTime/binStep); % bc times is given to us in seconds
                binRanges = 0:binStep:maxSpikeTime;
                
                % Set parameters for running assemply_patterns
                % Moved these 3 to outside of loop as they're constants
                %             opts.threshold.permutations_percentile = 95;
                %             opts.threshold.number_of_permutations = 20;
                %             opts.threshold.method = 'MarcenkoPastur';
                opts.Patterns.method = patternMethod;
                
                
                %%% Detect assemblies per ROI
                % also plot detected assesmblies and save those out
                
                for ROI = ROIs
                    
                    switch(char(ROI))
                        case 'ACA'
                            ROI_sessions = ACA_sessions;
                            ROI_color = [38, 237, 128]/255;
                        case 'CA1'
                            ROI_sessions = CA1_sessions;
                            ROI_color = [185, 5, 230]/255;
                        case 'CA3'
                            ROI_sessions = CA3_sessions;
                            ROI_color = [14, 141, 232]/255;
                    end
                    
                    [ROI_mat, ROI_AssemblyTemplates, n_assemblies_ROI] = detectAssemblies(ROI, ...
                        sessionIdx, sessionName, ROI_sessions, channels, ...
                        peakChannels, timesPerCluster, opts, binSize, ...
                        binRanges, binStep, base_path, patternMethod, ROI_color);
                end
                
                
            catch me
                keyboard
                % TO DO
                %     Activities = assembly_activity(AssemblyTemplates,ACA_mat);
                %
                %     figure(4),clf
                %     subplot(211)
                %     imagesc(ACA_mat)
                %     xlim([0 100])
                %     subplot(212)
                %     plot(Activities')
                %     xlim([0 100])
                
                badSessions = [badSessions, sessionIdx];
                
            end
        end
        
    end
    
end

%% Save out which sessions are bad to load in later if needed
% to not need to re-run the above for loop
% save('badSessions.mat', 'badSessions');
