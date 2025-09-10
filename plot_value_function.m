% plot_value_surface_3cases.m
% V(x1,0,x3,0) のサーフェスを 3 ケース比較（同一 Figure を再利用）

clear; clc;

% ===== 1) 比較する3ケース（ファイル名, 凡例） =====
cases = {
    'results_wo_codesign.mat', 'Without co-design';
    'results_SEPI.mat',       'System-Equivalence PI';
    'results_proposed.mat',    'Proposed co-design';
};

% ===== 2) 断面とグリッド設定（x2=0, x4=0） =====
x1_range = linspace(-1.5, 1.5, 41);
x3_range = linspace(-1.5, 1.5, 41);
[X1, X3] = meshgrid(x1_range, x3_range);
X2 = 0; X4 = 0;

% ===== 3) Figure を再利用 =====
fig = findobj('Type','figure','Tag','value_surface_compare_3');
if isempty(fig)
    fig = figure('Color','w','Position',[120 120 760 540], ...
                 'Tag','value_surface_compare_3');
else
    figure(fig); clf(fig);
end
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
view(ax, [-35 25]);  % 3D視点

% 色（お好みで変更可）
faceColors = [
    0.00 0.45 0.74;  % blue
    0.00 0.60 0.00;  % green
    0.85 0.33 0.10;  % orange/red
];

% ===== 4) 各ケースを計算・描画 =====
for i = 1:size(cases,1)
    file = cases{i,1}; leg = cases{i,2};
    V = eval_value_on_section(file, X1, X3, X2, X4);  % V(x1,0,x3,0)
    s = surf(ax, X1, X3, V, 'EdgeColor','none', ...
             'FaceAlpha', 0.75, 'FaceColor', faceColors(i,:));
    s.DisplayName = leg;
end

% ===== 5) 体裁 =====
title(ax, 'Comparison of the value functions');
xlabel(ax, 'x1 := x_L - y_d');
ylabel(ax, 'x3 := x_B');
zlabel(ax, 'V(x_1, 0, x_3, 0)');
legend(ax, 'Location','northeast');

out_png = 'value_function_comparison.png';
out_eps = 'value_fuction_comparison.eps';
exportgraphics(gcf, out_png, 'Resolution',300);
print(gcf, '-depsc', out_eps);


% ===== Helper: 断面上で V を評価 =====
function V = eval_value_on_section(matfile, X1, X3, X2, X4)
    S = load(matfile);

    x_sym   = ensure_sym(S.system.x);        % 1x4 sym
    Phi_sym = ensure_sym(S.setting.Phi);      % 1xN sym（必要なら cell→sym 変換）
    CN_last = S.outputs.CN_list(end,:).';     % N×1

    % V_N(x) = Phi(x) * C_N
    VN_sym = Phi_sym * CN_last;       % sym scalar

    % 関数化（スカラー引数）
    x1=x_sym(1); x2=x_sym(2); x3=x_sym(3); x4=x_sym(4);
    VN_fun = matlabFunction(VN_sym, 'Vars', {x1, x2, x3, x4});

    % グリッド上で評価（x2=0, x4=0）
    V = arrayfun(@(a,b) VN_fun(a, X2, b, X4), X1, X3);
end

function y = ensure_sym(y)
    if iscell(y), y = vertcat(y{:}); end
    if ~isa(y,'sym'), y = sym(y); end
end
