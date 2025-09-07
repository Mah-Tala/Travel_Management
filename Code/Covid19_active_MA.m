clc
clear all
close all

% -------------------------------------------------------------------------
% Active_MA.m
% Simulate active cases (asymptomatic + symptomatic) for B ∈ {15,20,22,25}.
% Loads:
%   ../../Datasets/Massachusetts_county/April/*.mat
%   ../../Datasets/Massachusetts_county/April/optimal_tau_B{15,20,22,25}.mat
% -------------------------------------------------------------------------

% Load base data
prefix = '../../Datasets/Massachusetts_county/April/';
d = dir([prefix,'*.mat']);
for i = 1:length(d)
    load([prefix, d(i).name]);
end
clear f_tau tau_optimal

tau0 = tau;
N    = population;
n    = length(N);

% Load optimized taus
load('../../Datasets/Massachusetts_county/April/optimal_tau_B15.mat'); tau_optimal_15 = tau_optimal; f_tau_15 = f_tau;
load('../../Datasets/Massachusetts_county/April/optimal_tau_B20.mat'); tau_optimal_20 = tau_optimal; f_tau_20 = f_tau;
load('../../Datasets/Massachusetts_county/April/optimal_tau_B22.mat'); tau_optimal_22 = tau_optimal; f_tau_22 = f_tau;
load('../../Datasets/Massachusetts_county/April/optimal_tau_B25.mat'); tau_optimal_25 = tau_optimal; f_tau_25 = f_tau;

% Plot f(tau)
figure;
plot(f_tau_15, 'LineWidth',2); hold on
plot(f_tau_20, 'LineWidth',2);
plot(f_tau_22, 'LineWidth',2);
plot(f_tau_25, 'LineWidth',2);
legend("B = 15","B = 20","B = 22","B = 25");
xlabel("# of iterations"); ylabel("f(\tau)"); grid on

% Build A for initial tau
A_init = tau0 * diag(1 ./ sum(diag(N) * tau0, 1)) * tau0' * diag(N);

% Params + beta calibration
s_init_all = sum(s_init)/n;
gamma_prime = 0.2; R_0 = 5.0;
gamma = gamma_prime; r_a = gamma_prime; r_s = gamma_prime;
epsilon = 0.32; alpha_transmission = 0.6754;

beta_hat = gamma_prime * R_0;
growth_rate_0 = s_init_all * beta_hat - gamma_prime;

beta_upper = 10; beta_low = 0;
while beta_upper > beta_low + 1e-5
    beta = (beta_upper + beta_low) / 2;
    M = [alpha_transmission * beta * diag(s_init) * A_init - (epsilon+r_a)* eye(n), beta * diag(s_init) * A_init;...
         epsilon * eye(n), -r_s * eye(n)];
    if max(real(eig(M))) - growth_rate_0 > 0
        beta_upper = beta;
    elseif max(real(eig(M))) - growth_rate_0 < 0
        beta_low = beta;
    else
        break
    end
end

% Build A for each optimized tau
A_optimal_15 = tau_optimal_15 * diag(1 ./ sum(diag(N) * tau_optimal_15, 1)) * tau_optimal_15' * diag(N);
A_optimal_20 = tau_optimal_20 * diag(1 ./ sum(diag(N) * tau_optimal_20, 1)) * tau_optimal_20' * diag(N);
A_optimal_22 = tau_optimal_22 * diag(1 ./ sum(diag(N) * tau_optimal_22, 1)) * tau_optimal_22' * diag(N);
A_optimal_25 = tau_optimal_25 * diag(1 ./ sum(diag(N) * tau_optimal_25, 1)) * tau_optimal_25' * diag(N);

parameter = [alpha_transmission beta epsilon r_a r_s];
tspan = 0:0.1:1000;

% ODE runs (3-compartment state: [s; x_a; x_s])
options = odeset('Events',@myevent);
x_prime_0 = [s_init; x_init; d_init];  % (kept as in your code)

[t_init, x_prime_init] = ode45(@(t1, x0) derivative_Covid19(t1, x0, A_init,       parameter), tspan, x_prime_0, options);
[t1,    x_prime_1   ]  = ode45(@(t1, x1) derivative_Covid19(t1, x1, A_optimal_15, parameter), tspan, x_prime_0, options);
[t2,    x_prime_2   ]  = ode45(@(t2, x2) derivative_Covid19(t2, x2, A_optimal_20, parameter), tspan, x_prime_0, options);
[t3,    x_prime_3   ]  = ode45(@(t3, x3) derivative_Covid19(t3, x3, A_optimal_22, parameter), tspan, x_prime_0, options);
[t4,    x_prime_4   ]  = ode45(@(t4, x4) derivative_Covid19(t4, x4, A_optimal_25, parameter), tspan, x_prime_0, options);

% Active = asymptomatic + symptomatic
alpha = 0.023;
x_active0 = sum(x_prime_0(n+1 : 3*n) .* [N ; N]) * exp(-alpha * tspan); % (kept as in your code)

x_active_init = sum(x_prime_init(:, n+1 : 3*n) .* repmat(N', [length(t_init), 2]), 2);
x_active1    = sum(x_prime_1(:,   n+1 : 3*n) .* repmat(N', [length(t1),    2]), 2);
x_active2    = sum(x_prime_2(:,   n+1 : 3*n) .* repmat(N', [length(t2),    2]), 2);
x_active3    = sum(x_prime_3(:,   n+1 : 3*n) .* repmat(N', [length(t3),    2]), 2);
x_active4    = sum(x_prime_4(:,   n+1 : 3*n) .* repmat(N', [length(t4),    2]), 2);

% Plot
figure
loglog(tspan, x_active0, 'Linestyle','-.', 'LineWidth', 1.5, 'Color', [0.49, 0.18, 0.56]); hold on
loglog(t4, x_active4, 'LineWidth', 1.5, 'Color', [0.31 0.31 0.31]);
loglog(t3, x_active3, 'LineWidth', 1.5, 'Color', [0.87 0.49 0   ]);
loglog(t2, x_active2, 'LineWidth', 1.5, 'Color', [0.47, 0.67, 0.19]);
loglog(t1, x_active1, 'LineWidth', 1.5, 'Color', [1 0 0]);
loglog(t_init, x_active_init, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
grid on
legend("$e ^ {-\alpha t}$", "B = 25", "B = 22", "B = 20", "B = 15", "initial $\tau (no lockdown)$","Interpreter","latex");
ylim([1e-1, 1e7]); xlim([1, 300]);
xlabel("Time (days)"); ylabel("Active cases");

% ----------------------------- local functions ---------------------------
function dxdt = derivative_Covid19(~, x_prime, A, parameter)
  [n,~] = size(A);
  alpha_transmission = parameter(1);
  beta     = parameter(2);
  epsilon  = parameter(3);
  r_a      = parameter(4);
  r_s      = parameter(5);
  s  = x_prime(1:n);
  % 3-block model as in your script ([s; x_a; x_s])
  Q = [zeros(3*n,n), ...
       [-alpha_transmission * beta * diag(s) * A ; ...
         alpha_transmission * beta * diag(s) * A - (epsilon + r_a)*eye(n); ...
         epsilon * eye(n)], ...
       [-beta * diag(s) * A; ...
         beta * diag(s) * A; ...
        -r_s * eye(n)]];
  dxdt = Q * x_prime;
end

function [value, isterminal, direction] = myevent(~, x_prime)
  value = min(x_prime); 
  isterminal = 1;  
  direction = 0;
end
