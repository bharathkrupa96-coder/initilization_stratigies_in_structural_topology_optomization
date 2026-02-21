% Create separate figures for each volume fraction
volfracs = [0.3 0.4 0.5 0.6 0.8];
strategy_names = {'uniform','random','quasi-sobol','structured','distance-constraint'};

for v = 1:length(volfracs)
    vf = volfracs(v);
    
    % Load data
    filename = sprintf('comparison_results%.1f.mat', vf);
    load(filename);
    
    % Create figure for this volume fraction
    figure('Name', sprintf('Optimized Structures - VF=%.1f', vf), ...
           'Position', [50 50 1400 800], 'NumberTitle', 'off');
    
    % Plot all 15 structures (5 strategies × 3 runs)
    for i = 1:15
        subplot(5, 3, i);
        
        % Plot the structure (black=material, white=void)
        imagesc(1 - multiresults(i).xPhys);
        colormap(gray);
        clim([0 1]);  % Ensure consistent color scale
        axis equal; 
        axis off;
        
        % Add title with strategy, run, and compliance
        title(sprintf('%s - Run %d\nC = %.3f', ...
            multiresults(i).strategy, ...
            multiresults(i).run, ...
            multiresults(i).c), ...
            'FontSize', 9, 'FontWeight', 'normal');
    end
    
    % Add overall title
    sgtitle(sprintf('Optimized Structures - Volume Fraction = %.1f', vf), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save the figure
    saveas(gcf, sprintf('Optimized_Structures_VF_%.1f.png', vf));
    
    fprintf('Created figure for VF = %.1f\n', vf);
end
