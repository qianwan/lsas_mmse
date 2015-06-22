function [x, Cx, xq, rh, eh] = iterative_cancellation_quantize_stats(L, M, K, H, y, N0, n, realx)
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
    x = zeros(L * K, 1);
    Cx = zeros(L * K, K);
    xq = zeros(L * K, 1);
    rh = [];
    eh = [];
    for p = 1 : L
        Hpp = H((p - 1) * M + 1 : p * M, (p - 1) * K + 1 : p * K);
        Up = U((p - 1) * M + 1 : p * M, :);
        yp = y((p - 1) * M + 1 : p * M) - sum(v((p - 1) * M + 1 : p * M, :), 2);
        x((p - 1) * K + 1 : p * K) = Hpp' / Up * yp;
        Cx((p - 1) * K + 1 : p * K, :) = eye(K) - Hpp' / Up * Hpp;
        for k = 1 : K
            s = x((p - 1) * K + k);
            re = real(s);
            im = imag(s);
            xq((p - 1) * K + k) = sign(re) / sqrt(2) + sign(im) / sqrt(2) * 1j;
            nom = exp(-abs(x((p - 1) * K + k) - xq((p - 1) * K + k))^2 / real(Cx((p - 1) * K + k, k)));
            den = 0;
            den = den + exp(-abs(x((p - 1) * K + k) - (+1/sqrt(2) + 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            den = den + exp(-abs(x((p - 1) * K + k) - (+1/sqrt(2) - 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            den = den + exp(-abs(x((p - 1) * K + k) - (-1/sqrt(2) - 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            den = den + exp(-abs(x((p - 1) * K + k) - (-1/sqrt(2) + 1j/sqrt(2)))^2 / real(Cx((p - 1) * K + k, k)));
            llr = log((nom / den) / (1 - nom / den));
            if realx((p - 1) * K + k) == xq((p - 1) * K + k)
                rh = [rh; llr];
            else
                eh = [eh; llr];
            end
        end
    end
