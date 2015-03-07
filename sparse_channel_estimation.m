function [h1, h2] = sparse_channel_estimation(y, S, M, N0, beta1, beta2, tau, T)
    h1 = beta1 * eye(M) / (tau * (beta1 + beta2) * eye(M) + N0 * eye(M)) * S' * y;
    h2 = beta2 * eye(M) / (tau * (beta1 + beta2) * eye(M) + N0 * eye(M)) * S' * y;
    h1 = h1 / norm(h1, 2) * sqrt(beta1 * M);
    h2 = h2 / norm(h2, 2) * sqrt(beta2 * M);
    SS = S' * S;
    for i = 1 : T
        mu = 2 * SS * h1 + 2 * SS * h2 - 2 * S' * y;
        h1 = h1 - 0.01 * mu;
        h2 = h2 - 0.01 * mu;
        h1 = h1 / norm(h1, 2) * sqrt(beta1 * M);
        h2 = h2 / norm(h2, 2) * sqrt(beta2 * M);
        yhat = S * (h1 + h2);
    end
