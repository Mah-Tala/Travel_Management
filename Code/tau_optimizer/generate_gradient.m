function [dlambda, lambda_0] = generate_gradient(tau, alpha_transmission, beta, s_init, r_a, r_s, epsilon, n, N, dz_numeric)
%GENERATE_GRADIENT  Gradient of lambda_max(M) wrt vec(tau)
% Inputs:
%   tau          : n×n current travel matrix
%   dz_numeric   : (2n*n^2)×(2n) stacked dM/dtau blocks
% Output:
%   dlambda      : n^2×1 gradient vector
%   lambda_0     : current dominant eigenvalue of M

    % current M at tau
    A_init = tau * diag(1 ./ sum(diag(N) * tau, 1)) * tau' * diag(N);
    Z0 = [alpha_transmission * beta * diag(s_init) * A_init - (epsilon + r_a) * eye(n),  beta * diag(s_init) * A_init; ...
          epsilon * eye(n),                                                             - r_s * eye(n)];

    % right eigenvector u_0 for dominant eigenvalue lambda_0
    [U, L] = eig(Z0);
    [~, idx] = max(diag(L));
    lambda_0 = L(idx, idx);
    u_0 = U(:, idx);

    % left eigenvector v_0 for the conjugate eigenvalue of Z0'
    Zs = Z0';
    [Evecs, Evals] = eig(Zs);
    lambda_bar = conj(lambda_0);
    [~, idx2] = min(abs(diag(Evals) - lambda_bar));
    v_0 = Evecs(:, idx2);
    v_star = v_0';   % row

    % block-diagonal V with v_star repeated n^2 times
    v_cell = mat2cell(repmat(v_star, 1, n^2), 1, ones(1, n^2) * (2 * n));
    V = blkdiag(v_cell{:});   % size: (n^2) × (2n*n^2)

    % gradient
    dlambda = (V * dz_numeric * u_0) / (v_star * u_0);   % n^2 × 1
end
