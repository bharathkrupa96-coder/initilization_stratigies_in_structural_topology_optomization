function ttest_all_volfrac_all_strategies_SAFE()

volfracs = [0.3 0.4 0.5 0.6 0.8];
strategy_names = {'uniform','random','quasi','structured','distance'};
alpha = 0.05;

fprintf('\n================ SAFE STATISTICAL ANALYSIS ================\n');

for v = 1:length(volfracs)
    vf = volfracs(v);
    filename = sprintf('comparison_results%.1f.mat', vf);
    load(filename);   % loads multiresults
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('Volume Fraction = %.1f\n', vf);
    fprintf('--------------------------------------------------\n');
    
    % Extract compliance
    comp = struct();
    for s = 1:length(strategy_names)
        strat = strategy_names{s};
        comp.(strat) = [multiresults(strcmp({multiresults.strategy}, strat)).c];%selects the structs whose strategy equals the current strat.
                                                                                % and get c of it 
    end
    
    % Pairwise comparison i stops at length - 1 so that j = i+1
    for i = 1:length(strategy_names)-1
        for j = i+1:length(strategy_names)
            
            s1 = strategy_names{i};
            s2 = strategy_names{j};
            
            c1 = comp.(s1);
            c2 = comp.(s2);
            
            % Remove NaNs
            c1 = c1(~isnan(c1));%gets only non nan values
            c2 = c2(~isnan(c2));
            
            % size and  Variance check
            if numel(c1) < 2 || numel(c2) < 2 || std(c1)==0 || std(c2)==0% mean -ind val ,square,sumof allsquares,res/n-1,final squrt(res)
                fprintf('%-12s vs %-12s | NOT TESTED (zero variance / deterministic)\n', ...
                    s1, s2);
                continue
            end
            
            % Valid t-test
            [h,p,~,stats] = ttest2(c1, c2, 'Alpha', alpha);
            
            fprintf('%-12s vs %-12s | p = %.4f | t = %+6.2f | %s\n', ...
                s1, s2, p, stats.tstat, ...
                ternary(h==1,'SIGNIFICANT','not significant'));
        end
    end
end

fprintf('\n================ END OF ANALYSIS =================\n');

end

function out = ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end

