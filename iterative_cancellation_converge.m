function [x, simu, calc] = iterative_cancellation_converge(L, M, K, H, y, N0, realx, n)
    v = zeros(L * M, L);
    W = zeros(L * M, L * M);
    W0 = zeros(L * M, L * M);
    U = zeros(L * M, M);
    simu = zeros(L * (n + 1), 1);
    calc = zeros(L * (n + 1), 1);
    for q = 1 : L
        %%%
        Hqq = H((q - 1) * M + 1 : q * M, (q - 1) * M + 1 : q * M);
        realxq = realx((q - 1) * K + 1 : q * K);
        yqq = y((q - 1) * M + 1 : q * M) - Hqq * realxq;
        calc(q) = norm(yqq, 2)^2 / M;
        %%%
        for p = 1 : L
            Hqp = H((q - 1) * M + 1 : q * M, (p - 1) * K + 1 : p * K);
            realxp = realx((p - 1) * K + 1 : p * K);
            W0((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) = Hqp * Hqp';
            %%%
            if p ~= q
                simu(q) = simu(q) + real(trace(W0((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M))) / M;
            end
            %%%
        end
    end
    for p = 1 : L
        Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        tmp = zeros(M);
        for l = 1 : L
            if l == p
                continue
            end
            tmp = tmp + W0((p - 1) * M + 1 : p * M, (l - 1) * M + 1 : l * M);
        end
        U((p - 1) * M + 1 : p * M, :) = Hpp * Hpp' + tmp + N0 * eye(M);
    end
    for t = 1 : n
        for p = 1 : L
            yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
            Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
            Up = U((p - 1) * M + 1: p * M, :);
            invUp = inv(Up);
            for q = 1 : L
                if q == p
                    continue
                end
                Hqp = H((q - 1) * M + 1 : q * M, (p - 1) * K + 1 : p * K);
                v((q - 1) * M + 1 : q * M, p) = Hqp * Hpp' * invUp * yp;
                W((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) = W0((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) - Hqp * Hpp' * invUp * Hpp * Hqp';
                simu(t * L + q) = simu(t * L + q) + real(trace(W((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M))) / M;
            end
        end
        for p = 1 : L
            Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
            Up = Hpp * Hpp' + N0 * eye(M);
            %%%
            yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
            realxp = realx((p - 1) * K + 1 : p * K);
            ypp = yp - Hpp * realxp;
            calc(t * L + p) = norm(ypp, 2)^2 / M;
            %%%
            for l = 1 : L
                if l == p
                    continue
                end
                Up = Up + W((p - 1) * M + 1 : p * M, (l - 1) * M + 1 : l * M);
            end
            U((p - 1) * M + 1: p * M, :) = Up;
        end
    end
    simu = simu + N0;
    x = zeros(L * K, 1);
    for p = 1 : L
        Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        Up = U((p - 1) * M + 1 : p * M, :);
        yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
        x((p - 1) * K + 1 : p * K) = Hpp' / Up * yp;
    end
