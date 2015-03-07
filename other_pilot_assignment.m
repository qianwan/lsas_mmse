function [A, B, C] = other_pilot_assignment(L, M, K, R, N0)
    A = zeros(L * K, 1);
    A(1 : K) = 1 : K;
    B = zeros(L * K, 1);
    B(1 : K) = 1 : K;
    for pilot = 1 : K
        for q = 2 : L
            S = 1 : q - 1;
            U = A((S - 1) * K + pilot)';
            SU = [S; U];
            minc = Inf;
            minu = NaN;
            for k = 1 : K
                if B((q - 1) * K + k) ~= 0
                    continue;
                end
                SUn = [SU, [q; k]];
                c = 0;
                for n = 1 : length(SUn)
                    b0 = SUn(1, n); u0 = SUn(2, n);
                    R0 = R((b0 - 1) * M + 1 : b0 * M, (b0 - 1) * K * M + (u0 - 1) * M + 1 : (b0 - 1) * K * M + u0 * M);
                    tmp = R0 + N0 / K * eye(M) / 3;
                    for m = 1 : length(SUn)
                        if n == m
                            continue;
                        end
                        b1 = SUn(1, m); u1 = SUn(2, m);
                        tmp = tmp + R((b0 - 1) * M + 1 : b0 * M, (b1 - 1) * K * M + (u1 - 1) * M + 1 : (b1 - 1) * K * M + u1 * M);
                    end
                    c = c + trace(R0 - R0^2 / tmp);
                end
                if c < minc
                    minc = c;
                    minu = k;
                end
            end
            B((q - 1) * K + minu) = pilot;
            A((q - 1) * K + pilot) = minu;
        end
    end
    C = 0;
    for pilot = 1 : K
        S = 1 : L;
        U = A((S - 1) * K + pilot)';
        SU = [S; U];
        c = 0;
        for n = 1 : length(SU)
            b0 = SU(1, n); u0 = SU(2, n);
            R0 = R((b0 - 1) * M + 1 : b0 * M, (b0 - 1) * K * M + (u0 - 1) * M + 1 : (b0 - 1) * K * M + u0 * M);
            tmp = R0 + N0 / K * eye(M);
            for m = 1 : length(SU)
                if n == m
                    continue;
                end
                b1 = SU(1, m); u1 = SU(2, m);
                tmp = tmp + R((b0 - 1) * M + 1 : b0 * M, (b1 - 1) * K * M + (u1 - 1) * M + 1 : (b1 - 1) * K * M + u1 * M);
            end
            c = c + trace(R0 - R0^2 / tmp);
        end
        C = C + c;
    end
    C = C / (L * M * K);
