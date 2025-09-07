function tau_min = quadratic_solver(tau_0, y_new, B, n)
%QUADRATIC_SOLVER  Project y_new onto { ||tau - tau_0||_1 <= B,  tau >= 0 }.
% Uses quadprog (Optimization Toolbox).

    % Decision var: [tau; z] ∈ R^(2n), with z ≥ 0 enforces L1 via linearization
    H = [eye(n), zeros(n); zeros(n), zeros(n)];
    f = -[y_new; zeros(n, 1)];

    % Inequalities A x <= b:
    % |tau - tau_0| <= z   →   tau - z <= tau_0,  -tau - z <= -tau_0
    A = zeros(2 * n, 2 * n);
    b = zeros(2 * n, 1);

    A(1:2:end, 1:n)      =  eye(n);   b(1:2:end, :) =  tau_0;
    A(2:2:end, 1:n)      = -eye(n);   b(2:2:end, :) = -tau_0;
    A(1:2:end, (n+1):2*n) = -eye(n);
    A(2:2:end, (n+1):2*n) = -eye(n);

    % Sum z_i <= B
    A = [A; zeros(1, n), ones(1, n)];
    b = [b; B];

    % Bounds: tau >= 0, z >= 0
    lb = zeros(2 * n, 1);

    % Solve
    x = quadprog(H, f, A, b, [], [], lb, [], []);
    tau_min = x(1:n);
end
