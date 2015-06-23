function [R, P] = generateReceiveCorrelation(L, M, K, BSs, UEs)
    R = zeros(L * M, L * K * M);
    P = 0;
    for q = 1 : L
        for p = 1 : L
            for k = 1 : K
                colOffset = (p - 1) * K * M + (k - 1) * M;
                d = abs(BSs(q) - UEs((p - 1) * K + k));
                intra_dist = abs(BSs(p) - UEs((p - 1) * K + k));
                tmp = 40 * log10(intra_dist / d);
                ps = 10^(tmp / 10);
                if p == q
                    P = P + ps * M;
                end
                R((q - 1) * M + 1 : q * M, colOffset + 1 : colOffset + M) = ps * eye(M);
            end
        end
    end
