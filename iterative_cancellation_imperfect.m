function x = iterative_cancellation_imperfect(L, M, K, Hhat, C, y, N0, n, detector)
    v = zeros(L * M, L);
    W = zeros(L * M, L * M);
    W0 = zeros(L * M, L * M);
    U = zeros(L * M, M);
    for q = 1 : L
        for p = 1 : L
            Hqp = Hhat((q - 1) * M + 1 : q * M, (p - 1) * K + 1 : p * K);
            W0((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) = Hqp * Hqp' + C((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M);
        end
    end
    for p = 1 : L
        Hpp = Hhat((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        tmp = zeros(M);
        for l = 1 : L
            if l == p
                continue
            end
            tmp = tmp + W0((p - 1) * M + 1 : p * M, (l - 1) * M + 1 : l * M);
        end
        if n ~= 0
            U((p - 1) * M + 1 : p * M, :) = Hpp * Hpp' + C((p - 1) * M + 1 : p * M, (p - 1) * M + 1 : p * M) + tmp + N0 * eye(M);
        else
            U((p - 1) * M + 1 : p * M, :) = Hpp * Hpp' + C((p - 1) * M + 1 : p * M, (p - 1) * M + 1 : p * M) + N0 * eye(M);
        end
    end
    for t = 1 : n
        for p = 1 : L
            yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
            Hpp = Hhat((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
            Up = U((p - 1) * M + 1: p * M, :);
            invUp = inv(Up);
            for q = 1 : L
                if q == p
                    continue
                end
                Hqp = Hhat((q - 1) * M + 1 : q * M, (p - 1) * K + 1 : p * K);
                v((q - 1) * M + 1 : q * M, p) = Hqp * Hpp' * invUp * yp;
                W((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) = W0((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) - Hqp * Hpp' * invUp * Hpp * Hqp';
            end
        end
        for p = 1 : L
            Hpp = Hhat((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
            Up = Hpp * Hpp' + C((p - 1) * M + 1 : p * M, (p - 1) * M + 1 : p * M) + N0 * eye(M);
            for l = 1 : L
                if l == p
                    continue
                end
                Up = Up + W((p - 1) * M + 1 : p * M, (l - 1) * M + 1 : l * M);
            end
            U((p - 1) * M + 1: p * M, :) = Up;
        end
    end
    x = zeros(L * K, 1);
    for p = 1 : L
        Hpp = Hhat((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        Up = U((p - 1) * M + 1 : p * M, :);
        yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
        if strcmp(detector, 'mmse') == 1 || strcmp(detector, 'MMSE') == 1
            x((p - 1) * K + 1 : p * K) = Hpp' / Up * yp;
        elseif strcmp(detector, 'zf') == 1 || strcmp(detector, 'ZF') == 1
            x((p - 1) * K + 1 : p * K) = pinv(Hpp) * yp;
        elseif strcmp(detector, 'mf') == 1 || strcmp(detector, 'MF') == 1
            x((p - 1) * K + 1 : p * K) = Hpp' * yp;
        end
    end
