function err = free_estimation_error(L, M, K, R, N0)
    err = 0;
    for q = 1 : L
        rowOffset = (q - 1) * M;
        for p = 1 : L
            if p ~= q
                continue
            end
            for k = 1 : K
                colOffset = (p - 1) * K * M + (k - 1) * M;
                r = R(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                t = trace(r - r^2 / (r + N0 * eye(M) / K / 3));
                err = err + t;
            end
        end
    end
    err = err / (L * M * K);
