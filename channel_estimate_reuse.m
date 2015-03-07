function [Hhat, C] = channel_estimate_reuse(L, M, K, H, R, A, N0, ratio)
    Hhat = H;
    C = zeros(L * M, L * M);
    K3 = K * ratio;
    for pilot = 1 : K3
        for q = 1 : L
            tR = zeros(1, L);
            sR = zeros(M, M);
            for p = 1 : L
                u = A((p - 1) * K + pilot);
                if u == 0
                    continue;
                end
                Rqpu  = R((q - 1) * M + 1 : q * M, (p - 1) * K * M + (u - 1) * M + 1 : (p - 1) * K * M + u * M);
                tR(p) = trace(Rqpu);
                sR    = sR + Rqpu;
            end
            [placeholder, index] = sort(tR, 'descend');
            toCan = zeros(M, M);
            for p = index
                u = A((p - 1) * K + pilot);
                if u == 0
                    continue;
                end
                Rqpu   = R((q - 1) * M + 1 : q * M, (p - 1) * K * M + (u - 1) * M + 1 : (p - 1) * K * M + u * M);
                Cqpu   = Rqpu - Rqpu^2 / (sR + toCan + N0 / K3 * eye(M));
                toCan  = toCan - Rqpu + Cqpu;
                tmp    = C((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M);
                C((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) = tmp + Cqpu;
                tmp = H((q - 1) * M + 1 : q * M, (p - 1) * K + u);
                Hhat((q - 1) * M + 1 : q * M, (p - 1) * K + u) = tmp + sqrt(Cqpu) * (randn(M, 1) + randn(M, 1) * 1j) / sqrt(2);
            end
        end
    end
