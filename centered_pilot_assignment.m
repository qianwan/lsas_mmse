function [A, C] = centered_pilot_assignment(L, M, K, R, N0)
    C = 0;
    A = zeros((L - 1) * K, 1);
    A(1 : K) = 1 : K;
    S = 1;
    for q = 2 : L - 1
        costMat = zeros(K, K);
        for k = 1 : K
            for pilot = 1 : K
                U = A((S - 1) * K + pilot)';
                SU = [S; U];
                SUn = [SU, [q; k]];
                c = 0;
                for n = 1 : length(SUn)
                    b0 = SUn(1, n); u0 = SUn(2, n);
                    R0 = R((L - 1) * M + 1 : L * M, (b0 - 1) * K * M + (u0 - 1) * M + 1 : (b0 - 1) * K * M + u0 * M);
                    tmp = R0 + N0 / K * eye(M) / 3;
                    for m = 1 : length(SUn)
                        if n == m
                            continue;
                        end
                        b1 = SUn(1, m); u1 = SUn(2, m);
                        tmp = tmp + R((L - 1) * M + 1 : L * M, (b1 - 1) * K * M + (u1 - 1) * M + 1 : (b1 - 1) * K * M + u1 * M);
                    end
                    c = c + trace(R0 - R0^2 / tmp);
                end
                costMat(k, pilot) = c;
            end
        end
        S = [S, q];
        [assignment, cost] = munkres(costMat);
        C = cost / (L - 1) / M / K;
        for pilot = 1 : length(assignment)
            A((q - 1) * K + pilot) = find(assignment(:, pilot));
        end
    end
