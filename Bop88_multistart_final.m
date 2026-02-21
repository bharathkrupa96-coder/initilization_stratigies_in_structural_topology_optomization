function [multiResults] = Bop88_multistart_final(volfrac)
% Bop88_multistart_final2 — Multi-start Top88 with 5 initialization strategies


%% ------------------- Input Handling -------------------
if nargin < 1
    volfrac = 0.5;
end

%% ------------------- User Settings ---------------------
nelx = 120;    
nely = 32;

penal = 3;
rmin = 1.5;
ft = 1;
tol = 1e-2;
maxiter = 200;
num_runs_per_strategy = 3;

seed = 0;
show_plots = true;

rng(seed);% reproduce

E0 = 1; Emin = 1e-9; nu = 0.3;

strategy_names = {'uniform','random','quasi-sobol','structured','distance-constraint'};
num_strategies = numel(strategy_names);

%% ------------------- Build FE Matrices -----------------% pre computed fea block
A11 = [12  3 -6 -3;   3 12  3  0;  -6  3 12 -3;  -3  0 -3 12];
A12 = [-6 -3  0  3;  -3 -6 -3 -6;   0 -3 -6  3;   3 -6  3 -6];
B11 = [-4  3 -2  9;   3 -4 -9  4;  -2 -9 -4 -3;   9  4 -3 -4];
B12 = [ 2 -3  4 -9;  -3  2  9 -2;   4  9  2  3;  -9 -2  3  2];
KE = 1/(1-nu^2)/24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11]);%sm of element 

nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
edofMat = zeros(nelx*nely, 8);%to store the global DOF numbers for every elemen


row = 1;
for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1) + ely;% lower left global node number
        n2 = (nely+1)*elx + ely;% lower right node
        edofMat(row,:) = [2*n1-1 2*n1 2*n2-1  2*n2 ...%xdoff ydoff of lower left right
                          2*(n2+1)-1 2*(n2+1) 2*(n1+1)-1 2*(n1+1)];% x, y doffs of upper,right,left
        row = row + 1;
    end
end

iK = reshape(kron(edofMat,ones(8,1))', 64*nelx*nely, 1);%row indecs
jK = reshape(kron(edofMat,ones(1,8))', 64*nelx*nely, 1);

%% ------------------- Boundary Conditions ---------------
F = sparse(2,1,-1, 2*(nely+1)*(nelx+1), 1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1), 2*(nelx+1)*(nely+1));
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

%% ------------------- Sensitivity Filter ----------------
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));

k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;% extracts id of current element
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;% extract id of neigbouring element in rmin
        k = k+1;
        iH(k) = e1; jH(k) = e2;
        sH(k) = max(0, rmin - sqrt((i1-i2)^2 + (j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);%% weights num of , row , num of col, sh corresponding
Hs = sum(H,2);%.

%% ------------------- Storage --------------------------- initilisation
totalRuns = num_strategies * num_runs_per_strategy;
template = struct('strategy','','run',0,'xPhys',[],...
                  'c',[],'iterations',[],'convergence',[],...
                  'vol_history',[],'runtime',[]);
multiResults = repmat(template,totalRuns,1);
all_final = cell(num_strategies, num_runs_per_strategy);

all_objectives = [];
all_runtimes = [];
all_mean_xPhys = [];
all_strategies = {};

%% ------------------- MULTI-START LOOP ------------------
idx = 0;
fprintf("\n================ MULTI-START OPTIMIZATION ================\n");
fprintf("Running %d strategies × %d runs each...\n\n", ...
         num_strategies, num_runs_per_strategy);

for s = 1:num_strategies
    stratName = strategy_names{s};

    fprintf("\n=============== STRATEGY %d: %s ===============\n", s, stratName);

    for run = 1:num_runs_per_strategy
        idx = idx + 1;
        fprintf("\n--- Starting Run %d/%d for Strategy %s ---\n", ...
                run, num_runs_per_strategy, stratName);

        %% -------------- INITIALIZATION -------------------
        switch s
            case 1
                x = volfrac * ones(nely,nelx);

            case 2
                x = volfrac * rand(nely,nelx);

            case 3
                try
                    p = sobolset(1,'Skip',1000,'Leap',200);%inbuilt matlab function
                    X = net(p, nelx*nely);% net returns the sobol points in p
                    x = volfrac * reshape(X, nely, nelx);
                catch
                    x = volfrac * rand(nely,nelx);
                end

            case 4
                x = volfrac * rand(nely,nelx);
                x(4:4:end,:) = 0;
                x(:,4:4:end) = 0;

            case 5
                fixedNodes = unique(ceil(fixeddofs/2));
                loadNode = ceil(find(F~=0)/2);
                critical = unique([fixedNodes(:); loadNode(:)]);

                pts = zeros(length(critical),2);
                for k = 1:length(critical)
                    nd = critical(k);
                    iy = mod(nd-1,nely+1) + 1;% extracts y cordinates
                    ix = floor((nd-1)/(nely+1)) + 1;% extract x cordinates
                    pts(k,:) = [ix iy];
                end

                [Xe,Ye] = meshgrid(1:nelx,1:nely);% creates 2 row matrix
                cx = Xe(:)+0.5; cy = Ye(:)+0.5;

                dmin = inf(length(cx),1);
                for k = 1:size(pts,1)
                    dmin = min(dmin, sqrt((cx-pts(k,1)).^2 + (cy-pts(k,2)).^2));% eularin distance formula used in filters
                end

                lambda = max(nelx,nely)/3;% length decay constant
                initMat = 0.05 + 0.95 * exp(-dmin/lambda);% exponential decay function 
                initMat = initMat * (volfrac / mean(initMat));% normalise to match volfrac
                initMat = reshape(initMat, nely, nelx);
                initMat = max(0.001, min(1, initMat));% to keep value in rance 0-1

                assignin('base', sprintf('initial_case5_matrix_run%d', run), initMat);
                x = initMat;
        end

        xPhys = x;

        %% ------------ OPTIMIZATION LOOP ------------------
        iter = 0; change = Inf;
        hist_c = []; hist_vol = [];
        tStart = tic;

        while change > tol && iter < maxiter
            iter = iter + 1;

            Ee = Emin + (xPhys(:)'.^penal)*(E0 - Emin);% simp equation
            sK = (KE(:)*ones(1,nelx*nely)) .* Ee;
            K = sparse(iK,jK,sK(:)); K = (K+K')/2;

            U(freedofs) = K(freedofs,freedofs)\F(freedofs);% equilibrium equation

            ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);%complience minimisation
            cval = sum(sum((Emin + xPhys.^penal*(E0-Emin)).*ce));%at element level

            dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;% derivative of sensitivity
            dv = ones(nely,nelx);% derivative of volume
         %% filters 
            if ft==1
                dc = reshape((H*(dc(:).*x(:)))./Hs ./ max(x(:),1e-9),nely,nelx);
            else
                dc = reshape(H*(dc(:)./Hs),nely,nelx);
                dv = reshape(H*(dv(:)./Hs),nely,nelx);
            end
%oc
            l1 = 0; l2 = 1e9; move = 0.2;

            while (l2-l1)/(l1+l2) > 1e-3
                lmid = 0.5*(l1+l2);
                xnew = max(0, max(x-move,min(1, min(x+move, x.*sqrt(-dc./dv/lmid)))));% update the design by dc,dv

                if ft==1
                    xPhys = xnew;
                else
                    xPhys(:) = (H*xnew(:))./Hs;
                end

                if mean(xPhys(:)) > volfrac
                    l1 = lmid;
                else
                    l2 = lmid;
                end
            end

            change = max(abs(xnew(:)-x(:)));
            x = xnew;

            hist_c(end+1) = cval;
            hist_vol(end+1) = mean(xPhys(:));

            if mod(iter,20)==0 || iter==1
                fprintf("  Iter %3d | Obj = %.4f | Vol = %.3f | Change = %.5f\n", ...
                        iter, cval, mean(xPhys(:)), change);
            end
        end

        runtime = toc(tStart);

        %% -------- Summary for this run --------
        fprintf(">> Completed Run %d | Strategy: %s | Obj = %.4f | Time = %.2f sec | Iter = %d\n", ...
                run, stratName, cval, runtime, iter);

        %% Save results
        result.strategy = stratName;
        result.run = run;
        result.xPhys = xPhys;
        result.c = cval;
        result.iterations = iter;
        result.convergence = hist_c;
        result.vol_history = hist_vol;
        result.runtime = runtime;

        multiResults(idx) = result;
        all_final{s,run} = xPhys;

        all_objectives(end+1) = cval;
        all_runtimes(end+1) = runtime;
        all_mean_xPhys(end+1) = mean(xPhys(:));
        all_strategies{end+1} = sprintf('%s-%d', stratName, run);

    end

    fprintf("\n===== End of Strategy %s =====\n", stratName);
end

fprintf("\n================= ALL RUNS COMPLETE =================\n");


%% ----------- Convergence History Plot -----------
if show_plots
    colors = lines(num_strategies);
    figure('Name','Convergence History','Units','normalized','Position',[0.1 0.1 0.8 0.7]);
    hold on;
    for idx = 1:length(multiResults)
        result = multiResults(idx);
        color_idx = strcmp(strategy_names, result.strategy);
        plot(result.convergence, 'Color', colors(color_idx,:), ...
             'LineWidth', 1.5, 'DisplayName', sprintf('%s-%d', result.strategy, result.run));
    end
    xlabel('Iteration');
    ylabel('Objective Value');
    title('Convergence History for All Runs');
    legend('Location','best', 'NumColumns', 3);
    grid on;
    xlim([0 maxiter]);
    saveas(gcf, 'Convergence_History.png');
end

%% ----------- Summary Statistics -----------
fprintf('\n=== SUMMARY STATISTICS ===\n');
for s = 1:num_strategies
    strat_mask = contains({multiResults.strategy}, strategy_names{s});
    strat_results = multiResults(strat_mask);
    
    objectives = [strat_results.c];
    runtimes = [strat_results.runtime];
    iterations = [strat_results.iterations];
    final_densities = all_mean_xPhys(strat_mask);
    
    fprintf('\n%s Strategy:\n', strategy_names{s});
    fprintf('  Objectives: %.4f ± %.4f (min: %.4f, max: %.4f)\n', ...
            mean(objectives), std(objectives), min(objectives), max(objectives));
    fprintf('  Runtimes: %.2f ± %.2f sec\n', mean(runtimes), std(runtimes));
    fprintf('  Iterations: %.1f ± %.1f\n', mean(iterations), std(iterations));
    fprintf('  Final Mean Density: %.4f ± %.4f\n', mean(final_densities), std(final_densities));
end


%% Export to workspace
assignin('base','top88_multi_final_results',multiResults);
assignin('base','all_objectives',all_objectives);
assignin('base','all_runtimes',all_runtimes);
assignin('base','all_mean_xPhys',all_mean_xPhys);
assignin('base','all_strategies',all_strategies);

fprintf("\nAll runs complete. Results stored in workspace.\n");
fprintf("Convergence history and summary statistics generated.\n");
fprintf("You can now generate Pareto plots using:\n");
fprintf(">> plot_pareto_curves()\n");
%% ----------- Plot Final Grid -----------
if show_plots
    figure('Name','Final Topologies','Units','normalized','Position',[0.05 0.05 0.9 0.85]);%b%t%L%R in %
    for s=1:num_strategies
        for r=1:num_runs_per_strategy
            subplot(num_strategies,num_runs_per_strategy,(s-1)*num_runs_per_strategy+r);
            imagesc(1-all_final{s,r}); colormap(gray); axis equal off;
            title(sprintf('%s run %d\nc = %.4f', ...
            strategy_names{s}, r, multiResults((s-1)*num_runs_per_strategy + r).c), ...
            'Interpreter','none', 'FontSize',10);
        end
    end
    saveas(gcf, 'All_Final_Topologies.png');
end

   
            

end

