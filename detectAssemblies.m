function [ROI_mat, ROI_AssemblyTemplates, n_assemblies_ROI] = detectAssemblies(ROI, ...
    session, sessionName, ROI_sessions, channels, peakChannels, timesPerCluster, ...
    opts, binSize, binRanges, binStep, base_path, patternMethod, ROI_color)

if ismember(session, ROI_sessions)
    
    % number of channels of the probe corresponding to each area
    ROI_channel_idx = find(strcmp(cellstr(channels.allen_ontology), ROI));
    
    % number of neurons recorded from each area
    is_ROI = ismember(peakChannels, ROI_channel_idx);
    
    ROI_timesPerCluster = timesPerCluster(is_ROI);
    n_ROI_clusters = length(ROI_timesPerCluster);
    
    % bin into 10 ms bins
    ROI_mat = nan(length(ROI_timesPerCluster),binSize-1);
    
    % if binSize is not -1
    % session 1 throws error bc leftside is 1 entry longer than right
    % side i.e. binSize = length(binCounts) + 1
    
    
    
    for cluster = 1:n_ROI_clusters
        currentCell = ROI_timesPerCluster(cluster);
        
        % 2nd dim is binSize
        binCounts = histcounts([currentCell{:}],binRanges);
        
        ROI_mat(cluster,:) = binCounts;
        
    end
    
    [ROI_AssemblyTemplates, n_assemblies_ROI] = assembly_patterns(ROI_mat,opts);
    
    % Plot detected assemblies
    
    for assembly = 1:3 % n_ROI_clusters
        
        figTitle = [sessionName, newline, 'Assembly Patterns in ', ROI, newline, ...
            opts.Patterns.method, ' with ',  opts.threshold.method, ' Treshold', ...
            newline, 'Assembly ', num2str(assembly), ' out of ', ...
            num2str(n_assemblies_ROI), newline, 'Bin Step = ', num2str(binStep)];
        
        figure;
        set(gcf, 'Visible', 'off');
        h = stem(ROI_AssemblyTemplates(:,assembly));
        title([figTitle{:}], 'Interpreter', 'None');
        prettifyFig();
        h.LineWidth = 2;
        h.MarkerFaceColor = ROI_color;
        h.Color = ROI_color;
        set(gca,'linewidth',3);
        
        % make 0-line centered
        axes = gca;
        set(gca, 'YLim', [-max(abs(axes.YLim)), max(abs(axes.YLim))]);

        % Save out plots
        figPath = [base_path, '/figures/AssemblyPatterns/', ...
            patternMethod, '/', num2str(binStep), '/', ...
            sessionName, '_', ROI '_assembly_patterns_', ...
            opts.Patterns.method, '_', opts.threshold.method, '_', 'binStep_', ...
            num2str(binStep), '_', num2str(assembly), '_out_of_', ...
            num2str(n_assemblies_ROI), '_assemblies.png'];
        
        saveas(gcf, [figPath{:}]);
        
        close all;
        
    end
    
end

end