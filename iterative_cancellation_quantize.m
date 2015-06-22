function [xhat, Cx, xq] = iterative_cancellation_quantize(L, M, K, H, y, N0, n, thllr)
    v = zeros(L * M, L);
    W = zeros(L * M, L * M);
    W0 = zeros(L * M, L * M);
    U = zeros(L * M, M);
    for q = 1 : L
        for p = 1 : L
            Hqp = H((q - 1) * M + 1 : q * M, (p - 1) * K + 1 : p * K);
            W0((q - 1) * M + 1 : q * M, (p - 1) * M + 1 : p * M) = Hqp * Hqp';
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
            end
        end
        for p = 1 : L
            Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
            Up = Hpp * Hpp' + N0 * eye(M);
            for l = 1 : L
                if l == p
                    continue
                end
                Up = Up + W((p - 1) * M + 1 : p * M, (l - 1) * M + 1 : l * M);
            end
            U((p - 1) * M + 1: p * M, :) = Up;
        end
    end
    xpl = zeros(L * K, 1);
    Cx = zeros(L * K, K);
    xq = zeros(L * K, 1);
    for p = 1 : L
        Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        Up = U((p - 1) * M + 1 : p * M, :);
        yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
        xpl((p - 1) * K + 1 : p * K) = Hpp' / Up * yp;
        Cx((p - 1) * K + 1 : p * K, :) = eye(K) - Hpp' / Up * Hpp;
        for k = 1 : K
            s = xpl((p - 1) * K + k);
            re = real(s);
            im = imag(s);
            xq((p - 1) * K + k) = sign(re) / sqrt(2) + sign(im) / sqrt(2) * 1j;
            nom = exp(-abs(xpl((p - 1) * K + k) - xq((p - 1) * K + k))^2 / real(Cx((p - 1) * K + k, k)));
            den = 0;
            den = den + exp(-abs(xpl((p - 1) * K + k) - (+1/sqrt(2) + 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            den = den + exp(-abs(xpl((p - 1) * K + k) - (+1/sqrt(2) - 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            den = den + exp(-abs(xpl((p - 1) * K + k) - (-1/sqrt(2) - 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            den = den + exp(-abs(xpl((p - 1) * K + k) - (-1/sqrt(2) + 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            llr = log((nom / den) / (1 - nom / den));
            if llr > thllr
                xpl((p - 1) * K + k) = xq((p - 1) * K + k);
                Cx((p - 1) * K + 1 : p * K, k) = 0;
                Cx((p - 1) * K + k, :) = 0;
            end
        end
    end
    xhat = zeros(L * K, 1);
    for p = 1 : L
        yp = y((p - 1) * M + 1 : p * M);
        Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        Ubp = Hpp * Hpp' + N0 * eye(M);
        for q = 1 : L
            Hpq = H((p - 1) * M + 1 : p * M, (q - 1) * K + 1 : q * K);
            if q == p
                continue
            end
            yp = yp - Hpq * xpl((q - 1) * K + 1 : q * K);
            Ubp = Ubp + Hqp * Cx((q - 1) * K + 1 : q * K, :) * Hqp';
        end
        xhat((p - 1) * K + 1 : p * K) = Hpp' / Ubp * yp;
    end
