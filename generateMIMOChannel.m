function [H, R] = generateMIMOChannel(L, M, BSs, K, UEs, R)
    H = zeros(L * M, L * K);
    for q = 1 : L
        rowOffset = (q - 1) * M;
        for p = 1 : L
            for k = 1 : K
                d = abs(BSs(q) - UEs((p - 1) * K + k));
                colOffset = (p - 1) * K * M + (k - 1) * M;
                r = R(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                H(rowOffset + 1 : rowOffset + M, (p - 1) * K + k) = sqrt(r) * (randn(M, 1) + randn(M, 1) * j) / sqrt(2);
            end
        end
    end
    return
