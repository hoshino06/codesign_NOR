function res = fun_codesign_once(alpha, UPDATE_RHO, rho_init, verbose, es)
% This function performs a single run of a control co-design (CCD) 
% iteration using Galerkin approximations of the Hamilton-Jacobi-Bellman 
% (HJB) framework. The procedure jointly updates the physical design 
% parameters (rho) and the feedback control law to minimize a composite 
% cost function consisting of structural and control performance terms. 
%
% Usage:
%   res = fun_codesign_once(alpha, UPDATE_RHO, rho_init, verbose)
%
% Input arguments:
%   alpha      - Weight on the structural cost (Js)
%   UPDATE_RHO - Boolean flag: if true, update design parameters rho;
%                if false, only update control law
%   rho_init   - Initial design parameter vector [mL, mB, dL, kB, dB]
%   verbose    - Verbosity level (0 = silent, 1 = iteration logs, 
%                2 = logs + save + plots)
%   es         - Setting of early stopping; if true, stop early
%
% Output:
%   res - A struct containing:
%         .timestamp   - Time of execution
%         .system      - System definitions (x, f, g, initial control)
%         .problem     - Problem settings (alpha, beta, w)
%         .setting     - Basis info, iteration count, update rate
%         .outputs     - Iteration histories of CN, J, rho
%
if nargin == 0
    UPDATE_RHO = false;  
    alpha    = 0;
    rho_init = [2,  20, 15, 15, 0.5]; % Initialize rho with the provided initial value
    verbose  = true;  % Set verbose mode for detailed output
    es       = false;
end
if nargin < 5
    es       = false;
end

if verbose >= 1
    fprintf('=== [%s] Codesign start ===\n', ...
        char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
end 

%% System equations

x = sym('x',[1 4]); 
f = @(p)[ ... % p = [mL, mB, dL, kB, dB]
    x(2); 
    -(1/p(1)+1/p(2))*p(3)*x(2) + p(4)/p(2)*x(3) + (0.2*p(4)/p(2))*x(3)^3 + (p(5)/p(2))*x(4); 
    x(4); 
    p(3)/p(2)*x(2) - p(4)/p(2)*x(3) - 0.2*p(4)/p(2)*x(3)^3 - p(5)/p(2)*x(4)];
g = @(p)[ 0; (1/p(1)+1/p(2)); 0; -1/p(2) ];

% Initial control law
u_init = -x(1);

% Basis functions
Phi_reduced = [... % 2nd and 4th order  
         x(1)^2, x(1)*x(2), x(1)*x(3), x(1)*x(4), ...
         x(2)^2, x(2)*x(3), x(2)*x(4), ...
         x(3)^2, x(3)*x(4), ...
         x(4)^2, ...
         x(1)^4, x(1)^3*x(2), x(1)^2*x(2)^2, x(1)*x(2)^3, x(2)^4, ...
         x(3)^4, x(3)^3*x(4), x(3)^2*x(4)^2, x(3)*x(4)^3, x(4)^4
         ];

Phi_full = [ ...
    % 2nd order (10)
    x(1)^2, x(1)*x(2), x(1)*x(3), x(1)*x(4), ...
    x(2)^2, x(2)*x(3), x(2)*x(4), ...
    x(3)^2, x(3)*x(4), ...
    x(4)^2, ...
    ... % 4th power (4)
    x(1)^4, x(2)^4, x(3)^4, x(4)^4, ...
    ... % 3+1 type (12)
    x(1)^3*x(2), x(1)^3*x(3), x(1)^3*x(4), ...
    x(2)^3*x(1), x(2)^3*x(3), x(2)^3*x(4), ...
    x(3)^3*x(1), x(3)^3*x(2), x(3)^3*x(4), ...＿
    x(4)^3*x(1), x(4)^3*x(2), x(4)^3*x(3), ...
    ... % 2+2 type (6)
    x(1)^2*x(2)^2, x(1)^2*x(3)^2, x(1)^2*x(4)^2, ...
    x(2)^2*x(3)^2, x(2)^2*x(4)^2, x(3)^2*x(4)^2, ...
    ... % 2+1+1 type (12)
    x(1)^2*x(2)*x(3), x(1)^2*x(2)*x(4), x(1)^2*x(3)*x(4), ...
    x(2)^2*x(1)*x(3), x(2)^2*x(1)*x(4), x(2)^2*x(3)*x(4), ...
    x(3)^2*x(1)*x(2), x(3)^2*x(1)*x(4), x(3)^2*x(2)*x(4), ...
    x(4)^2*x(1)*x(2), x(4)^2*x(1)*x(3), x(4)^2*x(2)*x(3), ...
    ... % 1+1+1+1 type (1)
    x(1)*x(2)*x(3)*x(4), ...
];

Phi = Phi_full;
N = numel(Phi);

%% CCD Problem Setting 
Js = @(p) 100*p(5)^2;
beta  = 1; 
rho   = rho_init;
rho_min = [1,  15, 10, 10, 0.1];
rho_max = [3,  25, 20, 20, 1.0];
num_iter = 10;
if UPDATE_RHO 
    rate = 0.05; % initial update rate of rho 
else
    rate = 0;
end

domain = [-1.5 1.5; -1.5 1.5; -1.5 1.5; -2 2];  % integration domain
stab_opts = struct( ...
            'gridN',      5, ...              % number of grid points per state axis (5^4 total points)
            'negTol',     1e-6, ...           % require dV < -negTol
            'eigTol',     1e-6, ...           % require eigenvalues real parts < -eigTol
            'vTol',       1e-9, ...           % tolerance for positive definiteness of V
            'bandRatio',  0.05);              % boundary band ratio for V_L level set

% calculation of the vector W
w = 1/int(int(int(int(1, x(1), domain(1,1), domain(1,2)), ...
                         x(2), domain(2,1), domain(2,2)), ...
                         x(3), domain(3,1), domain(3,2)), ...
                         x(4), domain(4,1), domain(4,2));  % weighting function
W = zeros(N,1);
for n = 1:N
    phi_n = Phi(n);
    funW  =  w * (phi_n);
    W(n,1) = int(int(int(int(funW, ...
                  x(1), domain(1,1), domain(1,2)), ...
                  x(2), domain(2,1), domain(2,2)), ...
                  x(3), domain(3,1), domain(3,2)), ...
                  x(4), domain(4,1), domain(4,2));
end

% Storage
CN_list  = nan(num_iter+1, N);
rho_list = nan(num_iter+1, length(rho));
J_list   = nan(num_iter+1, 3);

% Initial Evaluation
u = u_init;
l = 1000*x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + u^2;
[C_N, A, B] = galerkin_approx(x,f(rho),g(rho),u,l,Phi);       
J = alpha*Js(rho) + beta*(W.')*C_N;

rho_list(1,:)=rho; CN_list(1,:)=C_N; 
J_list(1,:)=[J, Js(rho), (W.')*C_N];
disp(['  rho=[' char(join(string(rho), ', ')) '] , J=' num2str(J), ' Jc=' num2str((W.')*C_N)])


%% Optimization Iteration
for i=1:num_iter
    if verbose >= 1
        fprintf('--- [%s] Iteration: %d / %d ---\n', ...
            char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), i, num_iter);
    end

    if UPDATE_RHO
        % ====== Normal mode: update ρ ======
        % Sensititvity wrt p 
        dAdp = sensitivity(x, f, g, u, Phi, rho); 
        dFdp = zeros(1, numel(rho));
        for k = 1:numel(rho)
            dAk = dAdp(:,:,k);
            dFdp(k) = alpha*Js(rho) + beta*(W.')*(A\dAk)*(A\B);
        end

        % Update controller
        V_N = Phi*C_N;
        u   = -(1/2)* g(rho)'*jacobian(V_N, x).';
        l   = 1000*x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + u^2;
    
        % Update p using Armijo backtracking with stability checks
        c1  = 1e-4;
        tau = 0.5;
        max_back = 20;
        J_old = J;
        accepted = false;
        if rate < 0.01; rate=0.01; end
        for ls = 1:max_back
            rho_trial = max(rho_min, min(rho - rate*dFdp, rho_max));
            [C_N, A, B] = galerkin_approx(x,f(rho_trial),g(rho_trial),u,l,Phi);
           
            J_trial = alpha*Js(rho_trial) + beta*(W.')*C_N;
            rhs = J_old + c1 * (dFdp * (rho_trial - rho).');
    
            [is_stable, cert] = stability_checks(x, f(rho_trial), g(rho_trial), u, V_N, domain, stab_opts);
    
            if J_trial <= rhs && is_stable
                % Accept
                accepted = true;
                break
            else
                rate = tau*rate;
            end
        end 
        if ~is_stable
            warning('Stability check failed after %d backtracks due to %s. Proceeding with smallest step.', max_back, cert.reason);
        elseif ~accepted
            warning('Armijo: step not accepted after %d backtracks. Proceeding with smallest step.', max_back);
        end

    else
        % ====== Comparison mode: fix ρ ======
        % Update controller
        V_N = Phi*C_N;
        u   = -(1/2)* g(rho)'*jacobian(V_N, x).';
        l   = 1000*x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + u^2;
        % Evaluation
        [C_N, A, B] = galerkin_approx(x,f(rho),g(rho),u,l,Phi);
        rho_trial = rho;
        J_trial   = alpha*Js(rho_trial) + beta*(W.')*C_N;
        [is_stable, cert] = stability_checks(x, f(rho_trial), g(rho_trial), u, V_N, domain, stab_opts);
        if ~is_stable
            warning('Stability check failed due to %s. Proceeding with smallest step.', cert.reason);
        end
    end

    rho = rho_trial;
    J   = J_trial;
    rho_list(i+1,:) = rho;
    J_list(i+1,:)   = [J Js(rho) (W.')*C_N];
    CN_list(i+1,:)  = C_N;
    fprintf('  rho=[%s], J=%.2f, Js=%.2f, Jc=%.2f, rate=%.2e, stability_check=%d\n', ...
        strjoin(string(rho), ', '), J, Js(rho), (W.')*C_N, rate, is_stable);

    if es && abs(J-J_list(i,1)) < 0.5
        fprintf('*** Early stopping at iter %d ***\n', i);
        break
    end

end

res = struct( ...
    'timestamp', datetime('now','Format','yyyy-MM-dd HH:mm:ss'), ...
    'system',   struct('x', x, 'f', f, 'g', g, 'u_init', u_init), ...
    'problem',  struct('alpha', alpha, 'beta', beta, 'w', w), ...
    'setting',  struct('Phi', Phi, 'N', N, 'num_iter', num_iter, 'rate', rate), ...
    'outputs',  struct('CN_list', CN_list, 'J_list', J_list, 'rho_list', rho_list) ...
);

if verbose >= 2
    fprintf('=== [%s] Save and plot results ===\n', ...
        char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    
    save('results.mat', 'res')
    
    plot(0:num_iter, log10(J_list(:,1)'), '-x' ) 
    xlabel('Iteration')
    ylabel('Log_{10}(J)')
end 

end % end of main function


%% Function: 4D integral helper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = integral4D(fun,x)
   I = int(int(int(int(fun, x(1), -1.5, 1.5), x(2), -1.5, 1.5), x(3), -1.5, 1.5), x(4), -2, 2);
end


%% Function: Galerkin Approximation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C_N, A, B] = galerkin_approx(x, f, g, u, l, Phi)
    N = numel(Phi);    
    dPhi = jacobian(Phi,x);
    A = zeros(N);
    B = zeros(N,1);
    parfor n = 1:N
        % ---- A ----
        for j = 1:N
            funA = (dPhi(j,:)*(f+g*u))*Phi(n);        
            A(n,j) = integral4D(funA,x);
        end
        % ---- B ----
        funB = l*Phi(n);
        B(n,1) = integral4D(funB,x);
    end
    C_N = -A\B;
end

%% Function: Sensitivity of A wrt p
function dAdp = sensitivity(x, f, g, u, Phi, p)
    N = numel(Phi);    
    dPhi = jacobian(Phi,x);
    % diff f and g wrt p
    P  = sym('P', [1, numel(p)]);
    df = subs(jacobian(f(P),P), P, p); % f wrt rho
    dg = subs(jacobian(g(P),P), P, p); % g wrt rho
    % dAdp
    dAdp = zeros(N, N, numel(p)); 
    for k = 1:numel(p)
        parfor n=1:N
            for j=1:N
                funDA = dPhi(j,:)*(df(:,k)+dg(:,k)*u) * Phi(n);
                dAdp(n,j,k) = integral4D(funDA,x);
            end
        end
    end
end


%% Function: Stability checks %%%%%%%%%%%%%%%
function [ok, cert] = stability_checks(x, f_sym, g_sym, u_sym, VN_sym, domain, opts)
% x: sym row [1xn], f_sym: sym nx1, g_sym: sym nxm, u_sym: sym (scalar) or m×1
% VN_sym: sym scalar Phi*C_N
% domain: nx2 [min max] for each state
% opts: struct with fields gridN, negTol, eigTol, vTol, bandRatio

    if ~isfield(opts,'gridN');     opts.gridN=5; end
    if ~isfield(opts,'negTol');    opts.negTol=1e-6; end
    if ~isfield(opts,'eigTol');    opts.eigTol=1e-6; end
    if ~isfield(opts,'vTol');      opts.vTol=1e-9; end
    if ~isfield(opts,'bandRatio'); opts.bandRatio=0.05; end

    n = numel(x);

    % ---- (i) Local eigenvalue test ----
    Acl_sym = jacobian(f_sym + g_sym*u_sym, x);    % ∂(f+g u)/∂x
    Acl = double(subs(Acl_sym, x, zeros(1,n)));    % at origin
    cert.Acl = Acl;
    cert.Acl_eigs = eig(Acl);
    if any(real(cert.Acl_eigs) >= -opts.eigTol)    % Require strictly negative real parts.
        ok = false;
        cert.reason = 'eig_failed';
        return
    end

    % ---- (ii) Local Lyapunov certification (via converse Lyapunov theorem) 
    % Construct converse Lyapunov
    Q = eye(n);
    try
        P = lyap(Acl', Q);  % A'P + P A = -Q
    catch
        ok = false;
        cert.reason = 'lyap_failed';
        return
    end
    if ~all(eig(P) > 0)
        ok = false;
        cert.reason = 'P_not_posdef';
        return
    end
    cert.P = P;
    VL_sym  = x*P*x.';  
    dVL_sym = jacobian(VL_sym, x) * (f_sym + g_sym*u_sym);
    VL_fun  = matlabFunction(VL_sym,  'Vars',{x});
    dVL_fun = matlabFunction(dVL_sym, 'Vars',{x});

    % Evaluation on sample grids
    [Xgrid, ~] = make_grid(domain, opts.gridN);
    VL_vals  = arrayfun(@(k) VL_fun(Xgrid(k,:)), 1:size(Xgrid,1))';
    dVL_vals = arrayfun(@(k) dVL_fun(Xgrid(k,:)), 1:size(Xgrid,1))';
    norm2    = sum(Xgrid.^2, 2);                          % ||x||^2

    % Search for a feasible sublevel set { V_L <= c }:
    % Heuristic: start from 90% of max(V_L) on the grid; if infeasible, halve c.
    VL_max  = max(VL_vals);
    c = 0.9*VL_max;
    ok_region = false;
    for shrink=1:20
        % On a boundary band around V_L = c, enforce \dot V_L <= -negTol * ||x||^2.
        band = opts.bandRatio * c;
        onB  = abs(VL_vals - c) <= band;
        if ~any(onB)
            % If the boundary sample is too sparse, check all interior points V_L<=c (excluding 0).
            onB = VL_vals <= c & VL_vals > 0;
        end
        if any(dVL_vals(onB) >= -opts.negTol.*max(norm2(onB), eps))
            c = c * 0.5;  % shrink c and retry
            continue
        else
            ok_region = true;
            break
        end
    end
    if ~ok_region
        ok = false;
        cert.reason = 'lyap_region_failed';
        cert.c = c;
        return
    end
    cert.c = c;
    cert.min_dV_on_boundary = max(dVL_vals(abs(VL_vals - c) <= opts.bandRatio*c));

    % ---- (iii) Outside the local ROA (V_L>c),  require dV_N <= -negTol*||x||^2 -----
    dVN_sym  = jacobian(VN_sym, x) * (f_sym + g_sym*u_sym);   % time derivative of V_N
    dVN_fun  = matlabFunction(dVN_sym, 'Vars',{x});
    dVN_vals = arrayfun(@(k) dVN_fun(Xgrid(k,:)), 1:size(Xgrid,1))';

    outside = VL_vals > c & norm2 > 1e-12; % exclude the origin
    if any(dVN_vals(outside) >= -opts.negTol.*max(norm2(outside), eps))
        ok = false;
        cert.reason = 'dVN_not_negative_outside';
        cert.max_dVN_outside = max(dVN_vals(outside));
        cert.c = c;
        return
    end

    ok = true;
    cert.reason = 'ok';
end


function [X, grids] = make_grid(domain, gridN)
    n = size(domain,1);
    grids = cell(1,n);
    for i=1:n
        grids{i} = linspace(domain(i,1), domain(i,2), gridN);
    end
    [G{1:n}] = ndgrid(grids{:});
    X = zeros(numel(G{1}), n);
    for i=1:n
        X(:,i) = G{i}(:);
    end
end
