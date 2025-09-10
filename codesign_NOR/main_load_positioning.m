%% Main script: run three scenarios and save results
% Assumes fun_codesign_once.m is on the path

clear; clc;

fprintf('=== [%s] Program start ===\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));

% ===== Common settings =====
alpha    = 0;                       % structural cost weight
rho_init = [2, 20, 15, 15, 0.5];    % initial design params [mL mB dL kB dB]
verbose  = 1;                       % 0: silent, 1: logs, 2: logs+save+plot (function internal)

% ===== 1) Proposed Codesign (update rho + controller) =====
UPDATE_RHO = true;
fprintf('\n=== Run: Proposed Codesign ===\n');
res_proposed = fun_codesign_once(alpha, UPDATE_RHO, rho_init, verbose);
save('results_proposed.mat', '-struct', 'res_proposed');


% ===== 2) System-Equivalence PI (SE-PI)  =====
UPDATE_RHO = false;
fprintf('\n=== Run: System-Equivalence PI (rho fixed with SEPI value) ===\n');
rho_SEPI = [1.13,  15, 10, 11.25, 0.375]; 
res_sepi = fun_codesign_once(alpha, UPDATE_RHO, rho_SEPI, verbose);
save('results_SEPI.mat', '-struct', 'res_sepi');

% ===== 3) Without Codesign (initial rho fixed) =====
UPDATE_RHO = false;
fprintf('\n=== Run: Without Codesign (rho fixed with initial value) ===\n');
res_wo =  fun_codesign_once(alpha, UPDATE_RHO, rho_init, verbose);
save('results_wo_codesign.mat','-struct', 'res_wo')

fprintf('\nSaved results:\n');
