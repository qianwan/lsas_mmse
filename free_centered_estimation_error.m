function [err, main] = free_centered_estimation_error(L, M, K, R, N0)
    err = 0;
    main = 0;
    p = L;
    for q = 1 : L - 1
    	for k = 1 : K
    		r = R((L - 1) * M + 1 : L * M, (q - 1) * K * M + (k - 1) * M + 1 : (q - 1) * K * M + k * M);
    		main = main + trace(r);
    		err = err + trace(r - r^2 / (r + N0 * eye(M) / K / 3));
    	end
    end
    err = err / (L - 1) / M / K;
    main = main / (L - 1) / M / K;
