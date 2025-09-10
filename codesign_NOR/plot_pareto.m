%% plot_pareto_from_saved.m
% Load saved results_*pareto*.mat and plot only the Pareto points
% connected by a curve. No dominated points, no legend, no title.

clear; clc;

% ---- Load saved file (robust to the earlier typo "parato") ----
fname_candidates = {'results_pareto.mat','results_parato.mat'};
loaded = false;
for k = 1:numel(fname_candidates)
    if exist(fname_candidates{k}, 'file')
        S = load(fname_candidates{k});  % expects variables: T, res_all
        fprintf('Loaded: %s\n', fname_candidates{k});
        loaded = true;
        break;
    end
end
if ~loaded
    error('results_pareto.mat (or results_parato.mat) not found in current folder.');
end

% ---- Extract data ----
T = S.T;  % table with columns: alpha, Js, Jc, J, rho_end
Js = T.Js(:);
Jc = T.Jc(:);

% Remove non-finite rows (just in case)
valid = isfinite(Js) & isfinite(Jc);
Js = Js(valid);
Jc = Jc(valid);

% ---- Compute Pareto front indices (minimize Js and Jc) ----
idx_pf = local_paretofront([Js, Jc]);

% ---- Sort the Pareto set by Js for a clean connecting curve ----
Js_pf = Js(idx_pf);
Jc_pf = Jc(idx_pf);
[Js_pf_sorted, ord] = sort(Js_pf);
Jc_pf_sorted = Jc_pf(ord);

% ---- Plot (Pareto points + connecting line) ----
fig = figure('Color','w','Position',[100 100 760 520]);
hold on; box on; grid on;

scatter(Js_pf_sorted, Jc_pf_sorted, 70, 'filled');      % Pareto points
plot(Js_pf_sorted, Jc_pf_sorted, '-', 'LineWidth', 1.5); % Connecting curve

xlabel('J_s (structural cost)');
ylabel('J_c (control cost)');

saveas(fig, 'pareto_plot.png');
print(fig, 'pareto_plot.eps', '-depsc');

% -------- Helper (local) --------
function idx = local_paretofront(X)
% X: N-by-2 matrix; minimize both columns
% Returns logical index of non-dominated points.
    N = size(X,1);
    idx = true(N,1);
    for a = 1:N
        if ~idx(a), continue; end
        xa = X(a,:);
        dom = all(bsxfun(@le, X, xa), 2) & any(bsxfun(@lt, X, xa), 2);
        dom(a) = false; % ignore self
        if any(dom)
            idx(a) = false;
        end
    end
end
