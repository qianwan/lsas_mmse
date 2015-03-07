clear;
L = 2;
M = 70;
K = 50;
SNRdB = 4;
numCases = 200;

SNR = 10^(SNRdB / 10);
N0 = 1 / SNR;

BSs = zeros(L, 1);
r = 2;
if L == 1
    BSs = 0;
elseif L == 2
    BSs = [0, r * 1j];
elseif L == 3
    BSs = [0, r * 1j, r * cos(pi / 6) + r * sin(pi / 6) * 1j];
elseif L == 4
    BSs = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
    F = [1, 2, 3, 3];
elseif L == 7
    BSs = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                r * cos(pi / 6) + r * (sin(pi / 6) - 1) * 1j, ...
                -r * 1j, ...
                -r * cos(pi / 6) + r * (sin(pi / 6) - 1) * 1j];
    F = [1, 2, 3, 3, 2, 3, 2];
end

err = 0;

for ci = 1 : numCases
    UEs = brownian(L, K, BSs, r / sqrt(3));
    [R, P] = generateReceiveCorrelation(L, M, K, BSs, UEs);
    H = generateMIMOChannel(L, M, BSs, K, UEs, R);
    xr = (randi([0, 1], L * K, 1) * 2 - 1) / sqrt(2);
    xi = (randi([0, 1], L * K, 1) * 2 - 1) / sqrt(2);
    x = xr + xi * 1j;
    y = H * x + (randn(L * M, 1) + randn(L * M, 1) * 1j) / sqrt(2) * sqrt(N0);

    %xhat = pinv(H) * y;
    %xhat = H' / (H * H' + N0 * eye(M)) * y;
    %xhat = iterative_cancellation(L, M, K, H, y, N0, 0);

    [A, cost, main, cros] = pilot_assignment(L, M, K, R, N0);
    [Hhat, C] = channel_estimate(L, M, K, H, R, A, N0);
    xhat = iterative_cancellation_imperfect(L, M, K, Hhat, C, y, N0, 2);
    err = err + sum(real(x) .* real(xhat) < 0) + sum(imag(x) .* imag(xhat) < 0);
end

fprintf(2, 'BER is %e\n', err / numCases / L / K / 2);
