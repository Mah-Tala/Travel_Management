% clc
% clear all
close all

% -------------------------------------------------------------------------
% optimize_tau_B20.m
% Projected gradient descent to minimize f(tau) = lambda_max(M) subject to
%   ||tau - tau0||_1 <= B,  tau >= 0
% Loads: population.mat, initial_rate_04_01.mat, travel.mat
% Saves: optimal_tau_B< B >.mat  (default B=20 → optimal_tau_B20.mat)
% -------------------------------------------------------------------------

% ---------- Load inputs (keep paths as in your pipeline) ------------------
prefix = '../../Datasets/Massachusetts_county/April/';
d = dir([prefix, '*.mat']);
for i = 1:length(d)
    load([prefix, d(i).name]);  % loads: population, s_init, tau, ...
end

N   = population;      % population vector (n×1)
n   = length(N);
tau0 = tau;            % starting point

clear tau

% ---------- Calibrate beta (bisection) to match initial growth rate ------
s_init_all   = sum(s_init) / n;
gamma_prime  = 0.2;
R_0          = 5.0;
beta_hat     = gamma_prime * R_0;
growth_rate_0 = s_init_all * beta_hat - gamma_prime;

% model parameters
gamma   = gamma_prime; 
r_a     = gamma_prime;
r_s     = gamma_prime;
epsilon = 0.32;
alpha_transmission = 0.6754;

A_init = tau0 * diag(1 ./ sum(diag(N) * tau0, 1)) * tau0' * diag(N);

beta_upper = 10; beta_low = 0; it_bis = 0;
while beta_upper > beta_low + 1e-5
    it_bis = it_bis + 1;
    beta = (beta_upper + beta_low) / 2;
    M = [alpha_transmission * beta * diag(s_init) * A_init - (epsilon + r_a) * eye(n),  beta * diag(s_init) * A_init; ...
         epsilon * eye(n),                                                           - r_s * eye(n)];
    if max(real(eig(M))) - growth_rate_0 > 0
        beta_upper = beta;
    elseif max(real(eig(M))) - growth_rate_0 < 0
        beta_low = beta;
    else
        break
    end
end

% ---------- Symbolic dA/dC for all entries (kept as in your code) --------
C = sym('C', n);
A = C * diag(1 ./ sum(diag(N) * C, 1)) * C' * diag(N);

dA_all = sym(zeros(n^3, n));
k = 0;
for i = 1:n
    for j = 1:n
        dA = diff(A, C(i, j));                    % n×n
        dA_all(n * k + 1 : n * (k + 1), :) = dA;  % stack blocks
        k = k + 1;
    end
end

% ---------- PGD Hyperparameters ------------------------------------------
B = 20;                         % L1 budget (default = 20)
step = 0.01;                    % base step size
step_size = step * ones(n^2, 1);
step_size(1:n+1:n^2, 1) = 0.6;  % larger steps on diagonal entries

tau_init  = reshape(tau0', [n^2, 1]);   % vectorized tau0
tolerance = 1e-5;
max_iter  = 800;
f_tau     = zeros(max_iter, 1);

% ---------- PGD loop ------------------------------------------------------
for ittr = 1:max_iter
    % (optional schedule)
    if ittr == 70
        step_size = 0.01 * ones(n^2, 1);
    end

    % build dA numerically at current tau0 (matrix form)
    dA_numeric = double(subs(dA_all, C, tau0));

    % build dM/dtau blocks → dz_numeric (2n × 2n) per variable, stacked
    dz_numeric = zeros(2 * n^3, 2 * n);
    for kk = 0:(n^2 - 1)
        blk = dA_numeric(n * kk + 1 : n * (kk + 1), :);
        dz_numeric(2 * n * kk + 1 : 2 * n * (kk + 1), :) = ...
            [alpha_transmission * beta * diag(s_init) * blk,  beta * diag(s_init) * blk; ...
             zeros(n),                                       zeros(n)];
    end

    % gradient of lambda_max(M) w.r.t. vec(tau) and current f value
    [gradient, f] = generate_gradient(tau0, alpha_transmission, beta, s_init, r_a, r_s, epsilon, n, N, dz_numeric);
    f_tau(ittr, 1) = f;

    % take a step in vector space
    tau_vec = reshape(tau0', [n^2, 1]);
    y_new   = tau_vec - step_size .* gradient;

    % project onto { ||tau - tau_init||_1 <= B,  tau >= 0 }
    if (norm(y_new - tau_init, 1) <= B) && (min(y_new) >= 0)
        proj = y_new;
    else
        proj = quadratic_solver(tau_init, y_new, B, numel(tau_init));
    end

    tau_new = proj;

    ittr    % print iteration

    % convergence
    if norm(tau_vec - tau_new) < tolerance
        display("reduction in tau less than tolerance");
        break;
    end

    % back to matrix form for next iteration
    tau0 = reshape(tau_new, [n, n])';
end

tau_optimal = tau0;

% ---------- Diagnostics plot ---------------------------------------------
figure; plot(f_tau, 'LineWidth', 2); grid on
xlabel('# of iterations'); ylabel('f(\tau)');

% ---------- Save (defaults to optimal_tau_B20.mat when B=20) -------------
save_name = sprintf('../../Datasets/Massachusetts_county/April/optimal_tau_B%d.mat', B);
save(save_name, 'f_tau', 'tau_optimal');
