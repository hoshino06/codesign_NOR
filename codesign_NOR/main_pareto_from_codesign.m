%% main_pareto_from_codesign.m
% Sweep alpha to trace the (Js, Jc) Pareto curve using fun_codesign_once.
% - Calls fun_codesign_once(alpha, true, rho_init, verbose)
% - Collects the final Js and Jc from res.outputs.J_list(end,:)
% - Saves results and plots a Pareto scatter/curve
%
% Note:
%   In your current fun_codesign_once, Js(p) == 0 (dummy). 
%   If you keep it as 0, all points will collapse on a vertical line.
%   Replace Js() in fun_codesign_once with your actual structural cost
%   before running this for a meaningful Pareto front.

clear; clc;

% -------- User settings --------
rho_init  = [1, 15, 10, 10, 1];     % initial rho (same length as in fun_codesign_once)
verbose   = 1;                      % 0/1/2 (2 => save & plot inside the function)
UPDATE_RHO = true;                  % co-design on
% Choose alpha grid (log-spaced is typical for scalarization sweeps)
alphas = logspace(-1, 1, 5);       % e.g., 1e-3 ... 1e3

% -------- Storage --------
numA = numel(alphas);
Js_vals  = nan(numA,1);
Jc_vals  = nan(numA,1);
J_vals   = nan(numA,1);
rho_last = nan(numA, numel(rho_init));
res_all  = cell(numA,1); 

fprintf('=== Pareto sweep start (N=%d alphas) ===\n', numA);
t_start = tic;

% -------- Sweep alpha --------
for i = 1:numA
    alpha = alphas(i);
    fprintf('\n[%2d/%2d] alpha = %.4g\n', i, numA, alpha);

    % Run co-design once for this alpha
    early_stopping = true;
    res = fun_codesign_once(alpha, UPDATE_RHO, rho_init, verbose, early_stopping);

    % Pull the last iteration results
    J_list   = res.outputs.J_list;   % columns: [J, Js, Jc]
    last_idx = find(all(~isnan(J_list),2), 1, 'last');
    last     = J_list(last_idx, :);

    J_vals(i)  = last(1);
    Js_vals(i) = last(2);
    Jc_vals(i) = last(3);
    rho_last(i,:) = res.outputs.rho_list(last_idx,:);

    % (Optional) Warm-start next run with the previous optimal rho
    rho_init = rho_last(i,:);
    res_all{i} = res;
end

elapsed = toc(t_start);
fprintf('\n=== Pareto sweep done. Elapsed: %.2f sec ===\n', elapsed);

% -------- Save results --------
T = table(alphas(:), Js_vals(:), Jc_vals(:), J_vals(:), ...
          arrayfun(@(k) {rho_last(k,:)}, 1:numA).', ...
          'VariableNames', {'alpha','Js','Jc','J','rho_end'});

save('results_pareto.mat', 'T', 'res_all');


% -------- Compute (approx.) Pareto front (non-dominated set) --------
idx_pf = local_paretofront([Js_vals, Jc_vals]);   % logical index

% -------- Plot --------
figure('Color','w','Position',[100 100 760 520]);
hold on; box on; grid on;

% Non-dominated set
scatter(Js_vals(idx_pf), Jc_vals(idx_pf), 70, 'filled', 'DisplayName','Pareto (non-dominated)');
% Dominated set
scatter(Js_vals(~idx_pf), Jc_vals(~idx_pf), 35, 'o', 'DisplayName','Dominated');

% Draw a connecting curve (sort by Js for readability)
[~,ord] = sort(Js_vals);
plot(Js_vals(ord), Jc_vals(ord), '-', 'LineWidth', 1.2, 'DisplayName','Sweep order');

% Label a few representative alphas
lbl_id = unique(round(linspace(1,numA,min(numA,6)))); % up to 6 labels
for k = lbl_id
    text(Js_vals(k), Jc_vals(k), sprintf('  \\alpha=%.2g', alphas(k)), ...
        'FontSize', 10, 'Interpreter','tex');
end

xlabel('J_s (structural cost)'); 
ylabel('J_c (control cost proxy, W^T C_N)');
title('Pareto Sweep by \alpha (co-design)');
legend('Location','bestoutside');

% ---- Save figure ----
saveas(fig, 'pareto_plot.png');
saveas(fig, 'pareto_plot.eps');

% -------- Helper (local) --------
function idx = local_paretofront(X)
% X: N-by-2 (minimize both columns)
% Returns logical index of non-dominated points
    N = size(X,1);
    idx = true(N,1);
    for a = 1:N
        if ~idx(a), continue; end
        xa = X(a,:);
        % dominated if ∃b: X(b,:) <= xa & any strict <
        dom = all(bsxfun(@le, X, xa), 2) & any(bsxfun(@lt, X, xa), 2);
        dom(a) = false; % ignore self
        if any(dom)
            idx(a) = false;
        end
    end
end
