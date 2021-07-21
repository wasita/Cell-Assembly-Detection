function [ROI_mat, ROI_AssemblyTemplates, n_assemblies_ROI] = detectAssemblies(ROI, ...
    session, sessionName, ROI_sessions, channels, peakChannels, timesPerCluster, ...
    opts, binSize, binRanges, binStep, base_path, patternMethod, ROI_color)

% number of channels of the probe corresponding to each area
ROI_channel_idx = find(strcmp(cellstr(channels.allen_ontology), ROI));

% number of neurons recorded from each area
is_ROI = ismember(peakChannels, ROI_channel_idx);

ROI_timesPerCluster = timesPerCluster(is_ROI);
n_ROI_clusters = length(ROI_timesPerCluster); % should match sum(is_ROI);

% bin into binSize ms bins
ROI_mat = nan(length(ROI_timesPerCluster),binSize);

for cluster = 1:n_ROI_clusters
    currentCell = ROI_timesPerCluster(cluster);
    
    % 2nd dim is binSize
    binCounts = histcounts([currentCell{:}],binRanges);
    
    ROI_mat(cluster,:) = binCounts;
    
end

[ROI_AssemblyTemplates, n_assemblies_ROI] = assembly_patterns(ROI_mat,opts);




% Try k-means clustering on AssemblyTemplates for ROI:
k = 2;
rng(1); % For reproducibility

for assembly = 1:size(ROI_AssemblyTemplates,2)
    
    idx = kmeans(ROI_AssemblyTemplates(:,assembly),k);
    
    figure;
    set(gcf, 'Visible', 'off');
    plot(ROI_AssemblyTemplates(idx==1,assembly),ROI_AssemblyTemplates(idx==1,assembly),'r.','MarkerSize',12)
    hold on
    plot(ROI_AssemblyTemplates(idx==2,assembly),ROI_AssemblyTemplates(idx==2,assembly),'b.','MarkerSize',12)
    
    legend('Cluster 1','Cluster 2', 'Cluster 3', ...
        'Location','NW')
    title(['Cluster Assignments', newline, 'Assembly ', num2str(assembly), ' out of ', ...
        num2str(n_assemblies_ROI), newline, 'k=', num2str(k)]);
    
    prettifyFig()
    
    % Save out plots
    figPath = [base_path, '/figures/AssemblyPatterns/', ...
        patternMethod, '/', num2str(binStep), '/', ROI, '/', ...
        sessionName, '_', patternMethod, '_kmeans_', 'k_', num2str(k), ...
        '_cluster_', num2str(assembly), '_out_of_', ...
        num2str(n_assemblies_ROI), '.png'];
    
    saveas(gcf, [figPath{:}]);
    
end




Plot detected assemblies

for assembly = 1:n_assemblies_ROI
    
    figTitle = [sessionName, newline, 'Assembly Patterns in ', ROI, newline, ...
        opts.Patterns.method, ' with ',  opts.threshold.method, ' Treshold', ...
        newline, 'Assembly ', num2str(assembly), ' out of ', ...
        num2str(n_assemblies_ROI), newline, 'Bin Step = ', num2str(binStep)];
    
    figure;
    set(gcf, 'Visible', 'off');
    h = stem(ROI_AssemblyTemplates(idx==1,assembly));
    title([figTitle{:}], 'Interpreter', 'None');
    prettifyFig();
    h.LineWidth = 2;
    h.MarkerFaceColor = ROI_color;
    h.Color = ROI_color;
    set(gca,'linewidth',3);
    
    
    %         hold on
    figure;
    set(gcf, 'Visible', 'off');
    h = stem(ROI_AssemblyTemplates(idx==2,assembly));
    h.MarkerFaceColor = [185, 5, 230]/255;
    set(gca,'linewidth',3);
    prettifyFig();
    h.LineWidth = 2;
    
    % make 0-line centered
    axes = gca;
    %     set(gca, 'YLim', [-max(abs(axes.YLim)), max(abs(axes.YLim))]);
    
    set(gca, 'YLim', [-max(abs(axes.YLim)), max(abs(axes.YLim))]);
    
    % Save out plots
    figPath = [base_path, '/figures/AssemblyPatterns/', ...
        patternMethod, '/', num2str(binStep), '/', ROI, '/', ...
        sessionName, '_', ROI '_assembly_patterns_', ...
        opts.Patterns.method, '_', opts.threshold.method, '_', 'binStep_', ...
        num2str(binStep), '_', num2str(assembly), '_out_of_', ...
        num2str(n_assemblies_ROI), '_assemblies.png'];
    
    saveas(gcf, [figPath{:}]);
    
    close all;
    
end

end