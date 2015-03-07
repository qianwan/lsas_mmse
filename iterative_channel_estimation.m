numCases = 1;
SNRdB = 6;
M = 100;
tau = 10;

SNR = 10^(SNRdB / 10);
N0 = 1 / SNR;

s = randn(tau, 1);
s = s / norm(s) * sqrt(tau);
S = kron(s, eye(M));

beta1 = 1;
beta2 = 1;

mse1 = 0;
mse2 = 0;
h1p = 0;
h2p = 0;

for ci = 1 : numCases
    h1 = sqrt(beta1) * eye(M) * (randn(M, 1) + randn(M, 1) * 1j) / sqrt(2);
    h2 = sqrt(beta2) * eye(M) * (randn(M, 1) + randn(M, 1) * 1j) / sqrt(2);
    y = S * (h1 + h2) + (randn(tau * M, 1) + randn(tau * M, 1) * 1j) / sqrt(2) * sqrt(N0);
    %hath1 = R1 / (tau * (R1 + R2) + N0 * eye(M)) * S' * y;
    %hath2 = R2 / (tau * (R1 + R2) + N0 * eye(M)) * S' * y;
    [hath1, hath2] = iterative_estimate_channel(y, beta1 * eye(M), beta2 * eye(M), S, M, tau, N0);
    %[hath1, hath2] = sparse_channel_estimation(y, S, M, N0, 1, 0.1, tau, 30);
    mse1 = mse1 + norm(h1 - hath1, 2)^2 / M;
    mse2 = mse2 + norm(h2 - hath2, 2)^2 / M;
    h1p = h1p + norm(h1, 2)^2 / M;
    h2p = h2p + norm(h2, 2)^2 / M;
end

fprintf(2, 'H1p is %f, MSE1 is %f\nH2p is %f, MSE2 is %f\n', 10 * log10(h1p / numCases), 10 * log10(mse1 / numCases), ...
        10 * log10(h2p / numCases), 10 * log10(mse2 / numCases));
