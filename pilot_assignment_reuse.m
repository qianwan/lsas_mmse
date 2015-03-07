function A = pilot_assignment_reuse(L, M, K, R, N0, F, ratio)
    K3 = K * ratio;
    A = zeros(L * K3, 1);
    if L <= 3 || L == 7
        for q = 1 : L
            ps = F(q);
            A((q - 1) * K3 + 1 + ps * K : (q - 1) * K3 + K + ps * K) = 1 : K;
        end
        return
    end
    if L == 4
        for q = 1 : 3
            ps = F(q);
            A((q - 1) * K3 + 1 + ps * K : (q - 1) * K3 + K + ps * K) = 1 : K;
        end
        costMat = zeros(K, K);
        for k = 1 : K
            for pilot = K * 2 + 1 : K * 3
                u = A((3 - 1) * K3 + pilot);
                cost = 0;
                for q = 1 : 1
                    Ru = R((q - 1) * M + 1 : q * M, (3 - 1) * K * M + (u - 1) * M + 1 : (3 - 1) * K * M + u * M);
                    Rk = R((q - 1) * M + 1 : q * M, (4 - 1) * K * M + (k - 1) * M + 1 : (4 - 1) * K * M + k * M);
                    if trace(Ru) > trace(Rk)
                        tmp = Ru - Ru^2 / (Ru + Rk + N0 / K3 * eye(M));
                        cost = cost + trace(tmp) + trace(Rk - Rk^2 / (tmp + Rk + N0 / K3 * eye(M)));
                    else
                        tmp = Rk - Rk^2 / (Ru + Rk + N0 / K3 * eye(M));
                        cost = cost + trace(tmp) + trace(Ru - Ru^2 / (tmp + Ru + N0 / K3 * eye(M)));
                    end
                end
                costMat(k, pilot - K * 2) = cost;
            end
        end
        [assignment, cost] = munkres(costMat);
        for pilot = 1 : length(assignment)
            A((4 - 1) * K3 + K * 2 + pilot) = find(assignment(:, pilot));
        end
        return
    end
