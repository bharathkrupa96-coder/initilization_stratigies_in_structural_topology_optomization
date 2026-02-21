function plot_combined_pareto_curves()
% PLOT_COMBINED_PARETO_ALL_VOLFRACS
% Color by strategy, marker by volume fraction

clearvars -except multiresults

fileList = {
    'comparison_results0.3.mat'
    
    'comparison_results0.4.mat'
    'comparison_results0.5.mat'
    'comparison_results0.6.mat'
    'comparison_results0.8.mat'
};
volfrac_list = [0.3 0.4 0.5 0.6 0.8];
 % This should match the function setting

% Containers
all_objectives = [];  all_runtimes = [];  all_mean_xPhys = [];
all_strategies = {};  all_volfracs = [];

fprintf('Loading data from %d files...\n', length(fileList));

for i = 1:length(fileList)%loops through each file in 'filelist'
    filename = fileList{i};% saves ith file in filename comparison_results0.3.mat
    vf = volfrac_list(i);%saves ith volfrac from 'volfrac_list' in vf
    
    if ~isfile(filename)% (~ , not operator) if the file do not exist the condition is true  
        warning('Missing: %s', filename); 
        continue; 
    end
    
    % Load the data - trying different variable names
    try
        data = load(filename);% loads .mat file stored in 'filename',,, data contains all variables in .mat 
        
        % Check which variable exists in the file
        if isfield(data, 'multiresults')% check if multiresults exist in the data 
            results = data.multiresults;%if yes save in results  
       
        else
            % Try to get the first variable in the file
            vars = fieldnames(data);% fieldnames returns all var name stored in data
            if ~isempty(vars)% not empty 
                results = data.(vars{1});
            else
                warning('No data found in %s', filename);%if empty display this
                continue;
            end
        end
        
        % Extract data from results and store seperatly in particular
        n = length(results);
        fprintf('Loaded %s (volfrac %.1f, %d runs)\n', filename, vf, n);
        
        for j = 1:n% store the data from 'results' in the variable respectively
            all_objectives(end+1) = results(j).c;
            all_runtimes(end+1) = results(j).runtime;
            all_mean_xPhys(end+1) = mean(results(j).xPhys(:));
            all_volfracs(end+1) = vf;
            
            % Determine strategy name and run number
            strategy = results(j).strategy;
            run_num = results(j).run;
            all_strategies{end+1} = sprintf('%s-%d', strategy, run_num);
        end
        
    catch ME % if file could not load print the warning message
        warning('Error loading %s: %s', filename, ME.message);
        continue;
    end
end

if isempty(all_objectives)% ensure the data is stored in 'all_objective'
    error('No data loaded! Check your .mat files.');
end

fprintf('\nTotal data points loaded: %d\n', length(all_objectives));
fprintf('Volume fractions found: %s\n', mat2str(unique(all_volfracs)));

% Get unique strategies NAME (extract base strategy name without run number)
base_strategies = {};
for i = 1:length(all_strategies)
    % Extract base strategy name (remove -runNumber)
    parts = strsplit(all_strategies{i}, '-');% split based on '-' as all stratigies contain non unique names n run number
    base_strat = parts{1};
    if ~any(strcmp(base_strategies, base_strat))
        base_strategies{end+1} = base_strat;% this var just have the stratigies name only 
    end
end

unique_volfracs = sort(unique(all_volfracs));

% Color map for strategies
cmap = lines(length(base_strategies));

% Markers for volume fractions
volfrac_markers = {'o', 's', '^', 'd', 'v', 'p', 'h', '<', '>', '*'};
if length(unique_volfracs) > length(volfrac_markers)
    error('Too many volume fractions! Add more markers.');
end

% =========================================================================
% Plot 1: Runtime vs Objective (Color: Strategy, Marker: Volume Fraction)
% =========================================================================
fig1 = figure('Name','Runtime vs Compliance — Color=Strategy, Marker=Volfrac','Position',[100 100 1000 700]);
hold on; grid on; box on;

% Plot points: color by strategy, marker by volume fraction
for s = 1:length(base_strategies)
    strat_name = base_strategies{s};
    strat_color = cmap(s,:);% extracts 1 row of colourmap to use as color
    
    for v = 1:length(unique_volfracs)
        vf = unique_volfracs(v);
        marker = volfrac_markers{v};
        
        % Find points for this strategy and volume fraction
        mask = contains(all_strategies, strat_name, 'IgnoreCase',true) & ...
               (all_volfracs == vf);
        
        if any(mask)
            scatter(all_runtimes(mask), all_objectives(mask), 100, ...
                    strat_color, marker, ...
                    'MarkerFaceColor', strat_color, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 1, ...
                    'HandleVisibility', 'off');
        end
    end
end

% Create legend entries for strategies (color only)
h_strat = [];
for s = 1:length(base_strategies)
    h_strat(s) = scatter(NaN, NaN, 100, cmap(s,:), 'o', ...
                         'MarkerFaceColor', cmap(s,:), ...
                         'MarkerEdgeColor', 'k', ...
                         'LineWidth', 1, ...
                         'DisplayName', base_strategies{s});
end

% Create legend entries for volume fractions (marker only)
h_volfrac = [];
for v = 1:length(unique_volfracs)
    h_volfrac(v) = scatter(NaN, NaN, 100, [0.3 0.3 0.3], volfrac_markers{v}, ...
                           'MarkerFaceColor', [0.7 0.7 0.7], ...
                           'MarkerEdgeColor', 'k', ...
                           'LineWidth', 1, ...
                           'DisplayName', sprintf('Vf = %.1f', unique_volfracs(v)));
end

% Global Pareto front
[pf_obj1, pf_x1] = find_pareto_points(all_objectives, all_runtimes);
[pf_x1_s, idx] = sort(pf_x1);
h_pareto = plot(pf_x1_s, pf_obj1(idx), 'k-', 'LineWidth', 1, ...
                'DisplayName', 'Global Pareto Front');

% Combine legends
legend_handles = [h_strat, h_volfrac, h_pareto];
legend_labels = [base_strategies, ...
                 arrayfun(@(x) sprintf('Vf = %.1f', x), unique_volfracs, 'UniformOutput', false), ...
                 {'Global Pareto Front'}];

% Set plot properties
xlabel('Runtime (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Compliance', 'FontSize', 14, 'FontWeight', 'bold');
title('Runtime vs Compliance — Strategy (Color) × Volume Fraction (Marker)', ...
      'FontSize', 16, 'FontWeight', 'bold');
legend(legend_handles, legend_labels, 'Location', 'bestoutside', 'FontSize', 11);
set(gca, 'FontSize', 12);
grid on;

% =========================================================================
% Plot 2: Mean Density vs Objective (Color: Strategy, Marker: Volume Fraction)
% =========================================================================
fig2 = figure('Name','Mean Density vs Compliance — Color=Strategy, Marker=Volfrac','Position',[100 100 1000 700]);
hold on; grid on; box on;

% Plot points
for s = 1:length(base_strategies)
    strat_name = base_strategies{s};
    strat_color = cmap(s,:);
    
    for v = 1:length(unique_volfracs)
        vf = unique_volfracs(v);
        marker = volfrac_markers{v};
        
        mask = contains(all_strategies, strat_name, 'IgnoreCase',true) & ...
               (all_volfracs == vf);
        
        if any(mask)
            scatter(all_mean_xPhys(mask), all_objectives(mask), 100, ...
                    strat_color, marker, ...
                    'MarkerFaceColor', strat_color, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 1, ...
                    'HandleVisibility', 'off');
        end
    end
end

% Create legend entries for strategies
h_strat2 = [];
for s = 1:length(base_strategies)
    h_strat2(s) = scatter(NaN, NaN, 100, cmap(s,:), 'o', ...
                          'MarkerFaceColor', cmap(s,:), ...
                          'MarkerEdgeColor', 'k', ...
                          'LineWidth', 1, ...
                          'DisplayName', base_strategies{s});
end

% Create legend entries for volume fractions
h_volfrac2 = [];
for v = 1:length(unique_volfracs)
    h_volfrac2(v) = scatter(NaN, NaN, 100, [0.3 0.3 0.3], volfrac_markers{v}, ...
                            'MarkerFaceColor', [0.7 0.7 0.7], ...
                            'MarkerEdgeColor', 'k', ...
                            'LineWidth', 1, ...
                            'DisplayName', sprintf('Vf = %.1f', unique_volfracs(v)));
end

% Global Pareto front for mean density
[pf_obj2, pf_x2] = find_pareto_points(all_objectives, all_mean_xPhys);
[pf_x2_s, idx] = sort(pf_x2);
h_pareto2 = plot(pf_x2_s, pf_obj2(idx), 'k-', 'LineWidth', 1, ...
                 'DisplayName', 'Global Pareto Front');

% Combine legends
legend_handles2 = [h_strat2, h_volfrac2, h_pareto2];
legend_labels2 = [base_strategies, ...
                  arrayfun(@(x) sprintf('Vf = %.1f', x), unique_volfracs, 'UniformOutput', false), ...
                  {'Global Pareto Front'}];

% Set plot properties
xlabel('Mean Density (mean(x_{Phys}))', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Compliance', 'FontSize', 14, 'FontWeight', 'bold');
title('Mean Density vs Compliance — Strategy (Color) × Volume Fraction (Marker)', ...
      'FontSize', 16, 'FontWeight', 'bold');
legend(legend_handles2, legend_labels2, 'Location', 'bestoutside', 'FontSize', 11);
set(gca, 'FontSize', 12);
grid on;

% Set reasonable x-axis limits
xlim_volfrac = [min(all_mean_xPhys)-0.02, max(all_mean_xPhys)+0.02];
set(gca, 'XLim', xlim_volfrac);

% =========================================================================
% Plot 3: Summary Statistics Table
% =========================================================================
fprintf('\n=====================================================\n');
fprintf('SUMMARY STATISTICS\n');
fprintf('=====================================================\n');
fprintf('%-15s %-10s %-12s %-12s %-10s\n', ...
        'Strategy', 'VolFrac', 'Min Comp', 'Max Comp', 'Avg Comp');
fprintf('%-15s %-10s %-12s %-12s %-10s\n', ...
        '---------------', '----------', '------------', '------------', '----------');

for s = 1:length(base_strategies)
    for v = 1:length(unique_volfracs)
        mask = contains(all_strategies, base_strategies{s}, 'IgnoreCase',true) & ...
               (all_volfracs == unique_volfracs(v));
        
        if any(mask)
            comps = all_objectives(mask);
            fprintf('%-15s %-10.1f %-12.4f %-12.4f %-10.4f\n', ...
                    base_strategies{s}, unique_volfracs(v), ...
                    min(comps), max(comps), mean(comps));
        end
    end
end

fprintf('\nBest overall solution: Compliance = %.4f\n', min(all_objectives));
fprintf('Volume fraction of best: %.1f\n', all_volfracs(find(all_objectives == min(all_objectives), 1)));

% Save the figures using the figure handles
saveas(fig1, 'runtime_vs_compliance.png');
saveas(fig2, 'density_vs_compliance.png');

fprintf('\nPlots saved as:\n');
fprintf('  - runtime_vs_compliance.png\n');
fprintf('  - density_vs_compliance.png\n');
fprintf('\nAnalysis complete!\n');

end
% pareto helper
% =========================================================================
function [pareto_obj, pareto_x] = find_pareto_points(obj, x)
    n = numel(obj);
    dominated = false(n,1);
    for i = 1:n
        for j = 1:n
            if i ~= j && (obj(j) <= obj(i) && x(j) <= x(i)) && (obj(j) < obj(i) || x(j) < x(i))
                dominated(i) = true;
                break;
            end
        end
    end
    pareto_obj = obj(~dominated);
    pareto_x   = x(~dominated);
end
