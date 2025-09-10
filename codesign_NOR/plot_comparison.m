
% ---- 1) Set your files and labels here ----
cases = {
    'results_wo_codesign.mat', 'Without Codesign',               'gs-';  % green square solid
    'results_SEPI.mat',       'System-Equivalence PI',          'b^-';  % blue triangle solid
    'results_proposed.mat',   'Proposed Codesign',            'ro-';  % red circle solid
};

% ---- 2) Loader ----
function J = load_J_series(matfile)
% Expects the MAT file to contain a variable `res` with:
%   res.outputs.J_list  (numeric vector or matrix)

    % --- 1) File check ---
    if ~isfile(matfile)
        error('File not found: %s', matfile);
    end
    % --- 2) Load and basic shape checks ---
    S = load(matfile);
    if isfield(S, 'outputs') && isstruct(S.outputs) && isfield(S.outputs, 'J_list')
        J = S.outputs.J_list(:,1); % First colmun of J_list
    elseif isfield(S, 'J_list')
        J = S.J_list(:,1);
    else
        error(['No J_list found. Expected either S.outputs.J_list or ', ...
               'top-level J_list in %s.'], matfile);
    end
end

% ---- 3) Load all series ----
numC = size(cases,1);
Js = cell(numC,1);
maxLen = 0;
for i = 1:numC
    Js{i} = load_J_series(cases{i,1});
    maxLen = max(maxLen, numel(Js{i}))-1;
end

% ---- 4) Plot ----
figure('Color','w','Position',[100 100 680 420]); hold on; grid on;
for i = 1:numC
    y = Js{i};
    x = 0:numel(y)-1;
    yplot = log10(max(y, realmin)); % guard for non-positives
    plt = plot(x, yplot, cases{i,3}, 'LineWidth',1.2, 'MarkerSize',6);
    % Optional: enforce color/marker to match style more strictly
    switch i
        case 1, set(plt, 'Color',[0.85 0.1 0.1]); % nicer red
        case 2, set(plt, 'Color',[0.1 0.4 0.9]);  % nicer blue
        case 3, set(plt, 'Color',[0.1 0.6 0.2]);  % nicer green
    end
end
xlabel('Number of iterations');
ylabel('Log_{10}(J)');
%title('Proposed vs System-Equivalence PI vs Without Codesign');
legend(cases(:,2), 'Location','best');
xlim([0 maxLen]);
box on;

% ---- 5) Save figures ----
out_png = 'policy_iteration_comparison.png';
out_eps = 'policy_iteration_comparison.eps';
exportgraphics(gcf, out_png, 'Resolution',300);
print(gcf, '-depsc', out_eps);

fprintf('Saved: %s, %s\n', out_png, out_eps);

% ----- 6) Read final rho and last/first valid  =====
function rho = load_last_rho(matfile)
    if ~isfile(matfile)
        error('File not found: %s', matfile);
    end
    S = load(matfile);
    if isfield(S,'outputs') && isstruct(S.outputs) && isfield(S.outputs,'rho_list')
        R = S.outputs.rho_list;
    elseif isfield(S,'rho_list')
        R = S.rho_list;
    else
        error('No rho_list found in %s', matfile);
    end
    if isvector(R), rho = R(:).'; else, rho = R(end,:); end
end

function v = last_valid(x)
    idx = find(isfinite(x) & x>0, 1, 'last');
    if isempty(idx), error('No positive/finite J found'); end
    v = x(idx);
end

Jbest_wo   = last_valid(Js{1});   % Without Codesign 
Jbest_sepi = last_valid(Js{2});   % SE-PI
Jbest_prop = last_valid(Js{3});   % Proposed
X = 100*(1 - Jbest_prop/Jbest_wo);         % vs initial controller design
Y = 100*(1 - Jbest_prop/Jbest_sepi); % vs SE-PI

% ===== Final costs, reductions, and optimized parameters =====
fprintf('\n=== Cost (final-iteration) ===\n');
fprintf('Without Codesign (final) : %.6g\n', Jbest_wo);
fprintf('SE-PI            (final) : %.6g\n', Jbest_sepi);
fprintf('Proposed         (final) : %.6g\n', Jbest_prop);

fprintf('\n=== Reductions for manuscript ===\n');
fprintf('vs initial controller : %.1f %%\n', X);
fprintf('vs SE-PI              : %.1f %%\n', Y);

% ---- Optimized parameters (final) table ----
rho_names = {'m_L','m_B','d_L','k_B','d_B'};
rho_sepi  = load_last_rho(cases{2,1});
rho_prop  = load_last_rho(cases{3,1});

fprintf('\n=== Optimized design parameters (final iteration) ===\n');
fprintf('%-6s | %-12s | %-12s\n', 'Var','SE-PI','Proposed');
fprintf('-------------------------------------------\n');
K = min([numel(rho_names), numel(rho_sepi), numel(rho_prop)]);
for k = 1:K
    fprintf('%-6s | %-12.6g | %-12.6g\n', rho_names{k}, rho_sepi(k), rho_prop(k));
end