% Compare best structure from each volume fraction
volfracs = [0.3 0.4 0.5 0.6 0.8];

figure('Name', 'Best Structures Comparison', ...
       'Position', [100 100 1600 400]);

for v = 1:length(volfracs)
    vf = volfracs(v);
    
    % Load data
    filename = sprintf('comparison_results%.1f.mat', vf);
    load(filename);
    
    % Find the run with minimum compliance (best structure)
    compliances = [multiresults.c];
    [min_c, best_idx] = min(compliances);
    best = multiresults(best_idx);
    
    % Plot
    subplot(1, 5, v);
    imagesc(1 - best.xPhys);
    colormap(gray);
    axis equal; axis off;
    
    % Add text annotation
    title(sprintf('VF = %.1f\nBest: %s-Run%d\nC = %.3f', ...
        vf, best.strategy, best.run, min_c), ...
        'FontSize', 11, 'FontWeight', 'bold');
    
    % Optional: Add boundary outline
    hold on;
    contour(1 - best.xPhys, [0.5 0.5], 'r-', 'LineWidth', 1);
    hold off;
end

sgtitle('Comparison of Best Structures Across Volume Fractions', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Save
saveas(gcf, 'Best_Structures_Comparison.png');

