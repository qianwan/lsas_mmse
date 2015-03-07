clear;
M = 20;
K = 10;
SNRdB = 0;
numCases = 100;

SNR = 10^(SNRdB / 10);
N0 = 1 / SNR;

BSs = 0;
F = 0;
Fr = 50;

r = 1;

err = 0;
merr = 0;

for ci = 1 : numCases
    UEs = brownian(1, K, BSs, r / sqrt(3));
    [R, P] = generateReceiveCorrelation(1, M, K, BSs, UEs);
    H = generateMIMOChannel(1, M, BSs, K, UEs, R);
    Xr = (randi([0, 1], K, Fr) * 2 - 1) / sqrt(2);
    Xi = (randi([0, 1], K, Fr) * 2 - 1) / sqrt(2);
    X = Xr + Xi * 1j;
    Y = H * X + (randn(M, Fr) + randn(M, Fr) * 1j) / sqrt(2) * sqrt(N0);
    A = pilot_assignment_reuse(1, M, K, R, N0, F, 1);
    [Hhat, C] = channel_estimate_reuse(1, M, K, H, R, A, N0, 1);
    Xhat = iterative_channel_estimate_detection(M, K, Hhat, Y, C, Fr, N0, 1);
    err = err + sum(sum(real(X) .* real(Xhat) < 0)) + sum(sum(imag(X) .* imag(Xhat) < 0));
    Xhat = Hhat' / (Hhat * Hhat' + C) * Y;
    merr = merr + sum(sum(real(X) .* real(Xhat) < 0)) + sum(sum(imag(X) .* imag(Xhat) < 0));
end

fprintf(2, 'BER is %e, with %d errors\n', err / numCases / K / Fr / 2, err);
fprintf(2, 'BER is %e, with %d errors\n', merr / numCases / K / Fr / 2, merr);
