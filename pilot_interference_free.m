L = 2;
M = 30;
K = 30;
SNRdB = 0;
numCases = 100;

SNR = 10^(SNRdB / 10);
N0 = 1 / SNR;
N0 = 1;

for K = 20 : 10 : 100
M = K;

BSs = [];

r = 2;
if L == 1
    BSs = 0;
elseif L == 2
    BSs = [0, r * 1j];
elseif L == 3
    BSs = [0, r * 1j, r * cos(pi / 6) + r * sin(pi / 6) * 1j];
elseif L == 3.5
    BSs = [0, r * 1j, r * cos(pi / 6) + r * sin(pi / 6) * 1j, r * cos(pi / 6) / 3 + r * sin(pi / 6) * 1j];
    L = 4;
elseif L == 4
    BSs = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * j];
elseif L == 7
    BSs = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * 1j, ...
                r * cos(pi / 6) + r * (sin(pi / 6) - 1) * 1j, ...
                -r * cos(pi / 6) + r * (sin(pi / 6) - 1) * 1j];
elseif L == 19

end

ferr  = 0;
merr  = 0; mmain = 0; mcros = 0;
oerr  = 0;

fcerr = 0; fmain = 0;
cerr  = 0;

for ci = 1 : numCases
    UEs = brownian(L, K, BSs, r / sqrt(3));
    [R, P] = generateReceiveCorrelation(L, M, K, BSs, UEs);
    H = generateMIMOChannel(L, M, BSs, K, UEs, R);
    %ferr = ferr + free_estimation_error(L, M, K, R, N0);
    [A, C, main, cros] = pilot_assignment(L, M, K, R, N0);
    merr = merr + C; mmain = mmain + main; mcros = mcros + cros;
    %[A, B, C] = other_pilot_assignment(L, M, K, R, N0);
    %oerr = oerr + C;

    %[err, main] = free_centered_estimation_error(L, M, K, R, N0);
    %fcerr = fcerr + err; fmain = fmain + main;
    %[A, C] = centered_pilot_assignment(L, M, K, R, N0);
    %cerr = cerr + C;
end

fprintf(2, 'Free   estimation error is %fdB, K = %d\n', 10 * log10(ferr  / numCases), K);
fprintf(2, 'My     estimation error is %fdB, K = %d, main is %fdB, cross is %fdB\n', 10 * log10(merr  / numCases), K, 10 * log10(mmain / numCases), 10 * log10(mcros / numCases));
fprintf(2, 'Other  estimation error is %fdB, K = %d\n', 10 * log10(oerr  / numCases), K);
fprintf(2, 'Freec  estimation error is %fdB, K = %d, main is %fdB\n', 10 * log10(fcerr / numCases), K, 10 * log10(fmain / numCases));
fprintf(2, 'Center estimation error is %fdB, K = %d\n', 10 * log10(cerr  / numCases), K);


end
