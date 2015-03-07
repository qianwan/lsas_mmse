function [A, C, main, cros] = pilot_assignment(L, M, K, R, N0)
    A = zeros(L * K, 1);
    A(1 : K) = 1 : K;
    S = [1];
    C = 0;
    main = 0;
    for q = 2 : L
        costMat = zeros(K, K);
        for k = 1 : K
            for pilot = 1 : K
                U = A((S - 1) * K + pilot)';
                SU = [S; U];
                SUn = [SU, [q; k]];
                c = 0;
                for n = 1 : length(SUn)
                    b0 = SUn(1, n); u0 = SUn(2, n);
                    R0 = R((b0 - 1) * M + 1 : b0 * M, (b0 - 1) * K * M + (u0 - 1) * M + 1 : (b0 - 1) * K * M + u0 * M);
                    tmp = R0 + N0 / K * eye(M);
                    for m = 1 : length(SUn)
                        if n == m
                            continue;
                        end
                        b1 = SUn(1, m); u1 = SUn(2, m);
                        tmp = tmp + R((b0 - 1) * M + 1 : b0 * M, (b1 - 1) * K * M + (u1 - 1) * M + 1 : (b1 - 1) * K * M + u1 * M);
                    end
                    c = c + trace(R0 - R0^2 / tmp);
                end
                costMat(k, pilot) = c;
            end
        end
        S = [S, q];
        [assignment, cost] = munkres(costMat);
        C = cost / L / M / K;
        for pilot = 1 : length(assignment)
            A((q - 1) * K + pilot) = find(assignment(:, pilot));
        end
    end
    cros = 0;
    for q = 1
        for pilot = 1 : K
            u0 = A((q - 1) * K + pilot);
            tmp  = 0;
            tmp2 = 0;
            for p = 2 : L
                u1 = A((p - 1) * K + pilot);
                r  = R((p - 1) * M + 1 : p * M, (q - 1) * K * M + (u0 - 1) * M + 1 : (q - 1) * K * M + u0 * M);
                R1 = R((p - 1) * M + 1 : p * M, (p - 1) * K * M + (u1 - 1) * M + 1 : (p - 1) * K * M + u1 * M);
                tmp  = tmp  + trace(r);
                residual = R1 - R1^2 / (R1 + r + N0 / K * eye(M));
                tmp2 = tmp2 + trace(r - r^2 / (r + residual + N0 / K * eye(M)));
            end
            tmp  = tmp  / (L - 1);
            tmp2 = tmp2 / (L - 1);
            main = main + tmp;
            cros = cros + tmp2;
        end
    end
    main = main / M / K;
    cros = cros / M / K;
