antennas = [20, 40, 60, 80, 100];
free = [-78.6, -79.3, -79.9, -80.4, -80.9];
other2 = [-70.5, -72.3, -73.5, -74.3, -75.1] - 3.5;
other3 = [-62.1, -64.4, -66.2, -67.3, -68.1] - 6.7;
self2 = [-77.51, -78.59, -79.65, -80.31, -80.81] + 0.4;
self3 = [-76.23, -77.74, -78.98, -79.74, -80.42] + 0.4;
plot(antennas, free, '^-', antennas, other2, '-o', antennas, other3, '-s', antennas, self2, '-d', antennas, self3, '-v', 'LineWidth', 1.5);
legend('Interference free scenario', 'User section, L = 2', 'User section, L = 3', 'Cell-wise greedy assignment, L = 2', 'Cell-wise greedy assignment, L = 3');
xlabel('Number of antennas (K = M)');
ylabel('Estimation error (dB)');
axis([20, 100, -81, -60]);
grid on;


dB = [-10, -8, -6, -4, -2, 0, 2, 4];
% L = 1, M = K = 100, perfect channel estimation
mmselimit = [5.0e-2, 3.1e-2, 1.7e-2, 8.7e-3, 3.9e-3, 1.6e-3, 5.6e-4, 1.7e-4];
zflimit = [0.38, 0.359, 0.332, 0.298, 0.252, 0.206, 0.165, 0.139];
% L = 2, M = K = 100, imperfect channel estimation
cell2limit = [6.61e-2, 4.99e-2, 3.63e-2, 2.75e-2, 2.13e-2, 1.71e-2, 1.44e-2, 1.33e-2];
cell2limitIter1 = [5.71e-2, 3.73e-2, 2.30e-2, 1.38e-2, 7.58e-3, 4.38e-3, 2.44e-3, 1.41e-3];
cell2limitIter2 = [5.67e-2, 3.82e-2, 2.28e-2, 1.31e-2, 7.03e-3, 3.55e-3, 1.71e-3, 8.49e-4];
semilogy(dB, mmselimit, '-^', dB, zflimit, '-s', dB, cell2limit, '-d', dB, cell2limitIter1, '-o', dB, cell2limitIter2, '-v', 'LineWidth', 1.5);
legend('Single cell, MMSE, perfect channel estimation', 'Single cell, ZF, perfect channel estimation');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

dB = [-10, -8, -6, -4, -2, 0, 2, 4];
% L = 7, M = 100, K = 90, imperfect channel estimation
cell7limit = [8.98e-2, 7.47e-2, 6.41e-2, 5.86e-2, 5.32e-2, 5.11e-2, 4.75e-2, 4.62e-2];
cell7limitIter1 = [6.06e-2, 4.26e-2, ];
cell7limitIter2 = [];
semilogy(dB, cell7limit, '-^', [-10:2:4], cell7limitIter1, '-s', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('BER')
grid on;



% L = 7, M = 100
users = 76 : 4 : 100;
Iter0M4 = [3.09e-2, 3.90e-2, 4.52e-2, 5.43e-2, 6.33e-2, 7.06e-2, 7.84e-2];
Iter1M4 = [7.42e-3, 1.09e-2, 1.75e-2, 2.30e-2, 3.13e-2, 3.78e-2];
Iter2M4 = [4.29e-3, 8.41e-3, 1.21e-2, 1.77e-2, 2.50e-2, 3.20e-2];
Iter0M0 = [3.00e-2, 3.78e-2, 4.58e-2, 5.36e-2, 6.16e-2, 7.06e-2];
Iter0M1 = [1.97e-3, 4.05e-3, 8.06e-3, 1.29e-2, 1.90e-2, 2.66e-2];
Iter0M2 = [3.95e-4, 1.26e-3, 3.34e-3, 7.08e-3, 1.26e-2, 1.96e-2];
semilogy(users, Iter0M4, '-^', users, Iter1M4, '-s', users, Iter2M4, '-d', users, Iter0M0, '-o', users, Iter0M1, '-v', users, Iter0M2, '-*', 'LineWidth', 1.5);
legend('SNR = -4dB, #0', 'SNR = -4dB, #1');
xlabel('Number of Users (K)');
ylabel('BER');
grid on;









