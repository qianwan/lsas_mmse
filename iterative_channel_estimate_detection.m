function Xhat = iterative_channel_estimate_detection(M, K, Hhat, Y, C, Fr, N0, n)
    Xhat = Hhat' / (Hhat * Hhat' + C) * Y;
    if n == 0
        return;
    end
    Cd = zeros(M, M * K);
    for k = 1 : K
        Cd(:, (k - 1) * M + 1 : k * M) = C / K;
    end
    for ni = 1 : n
        Exx = eye(K) - Hhat' / (Hhat * Hhat' + C) * Hhat;
        dExx = diag(Exx);
        Yt = Y - Hhat * Xhat;
        yv = reshape(Yt, [], 1);
        Eyy = zeros(Fr * M, Fr * M);
        for k = 1 : K
            xhatk = Xhat(k, :);
            Eyy = Eyy + kron(xhatk.' * conj(xhatk), Cd(:, (k - 1) * M + 1 : k * M));
            for i = 1 : K
                Eyy = Eyy + kron(Exx(k, i) * eye(Fr), Hhat(:, k) * Hhat(:, i)');
            end
            Eyy = Eyy + kron(Exx(k, k) * eye(Fr), Cd(:, (k - 1) * M + 1 : k * M));
        end
        Eyy = Eyy + N0 * eye(Fr * M);
        Hhatnew = Hhat;
        Cnew = zeros(size(C));
        Cdnew = zeros(size(Cd));
        tmp = zeros(M, Fr * M);
        Hdelta = zeros(size(Hhat));
        for k = 1 : K
            if k > 1
                xhatk_1 = Xhat(k - 1, :);
                Xhatk_1 = kron(xhatk_1.', eye(M));
                yv = yv - Xhatk_1 * Hdelta(:, k - 1);
                Eyy = Eyy + Xhatk_1 * (Cdnew(:, (k - 2) * M + 1 : (k - 1) * M) - Cd(:, (k - 2) * M + 1 : (k - 1) * M)) * Xhatk_1';
            end
            xhatk = Xhat(k, :);
            Xhatk = kron(xhatk.', eye(M));
            Edhky = Cd(:, (k - 1) * M + 1 : k * M) * Xhatk';
            Hdelta(:, k) = Edhky / Eyy * yv;
            Hhatnew(:, k) = Hhatnew(:, k) + Hdelta(:, k);
            Cdnew(:, (k - 1) * M + 1 : k * M) = Cd(:, (k - 1) * M + 1 : k * M) - Edhky / Eyy * Edhky';
            Cnew = Cnew + Cdnew(:, (k - 1) * M + 1 : k * M);
        end
        Xhat = Hhatnew' / (Hhatnew * Hhatnew' + Cnew) * Y;
        if ni == n
            break;
        end
        Hhat = Hhatnew;
        C = Cnew;
        Cd = Cdnew;
    end
