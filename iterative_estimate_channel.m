function [h1, h2] = iterative_estimate_channel(y, R1, R2, S, M, tau, N0)
    h1 = R1 / (tau * (R1 + R2) + N0 * eye(M)) * S' * y;
    h2 = R2 / (tau * (R1 + R2) + N0 * eye(M)) * S' * y;
    C1 = R1 - R1^2 / (R1 + R2 + N0 * eye(M) / tau);
    C2 = R2 - R2^2 / (R1 + R2 + N0 * eye(M) / tau);
    h1p = R1 / (tau * (R1 + C2) + N0 * eye(M)) * S' * (y - S * h2);
    h2p = R2 / (tau * (R2 + C1) + N0 * eye(M)) * S' * (y - S * h1);
    C1p = R1 - R1^2 / (R1 + C2 + N0 * eye(M) / tau);
    C2p = R2 - R2^2 / (R2 + C1 + N0 * eye(M) / tau);
    h1 = h1p; C1 = C1p;
    h2 = h2p; C2 = C2p;

    h1p = R1 / (tau * (R1 + C2) + N0 * eye(M)) * S' * (y - S * h2);
    h2p = R2 / (tau * (R2 + C1) + N0 * eye(M)) * S' * (y - S * h1);
    C1p = R1 - R1^2 / (R1 + C2 + N0 * eye(M) / tau);
    C2p = R2 - R2^2 / (R2 + C1 + N0 * eye(M) / tau);
    h1 = h1p; C1 = C1p;
    h2 = h2p; C2 = C2p;    
