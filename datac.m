% Warning, obsolete
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



dB = -10 : 2 : 4;
% L = 1, M = K = 100, perfect channel estimation
mmselimit = [5.0e-2, 3.1e-2, 1.7e-2, 8.7e-3, 3.9e-3, 1.6e-3, 5.6e-4, 1.7e-4];
zflimit   = [0.38, 0.359, 0.332, 0.298, 0.252, 0.206, 0.165, 0.139];
% L = 2, M = K = 100, perfect channel estimation
%cell2Iter0Perf   = [6.61e-2, 4.99e-2, 3.63e-2, 2.75e-2, 2.13e-2, 1.71e-2, 1.44e-2, 1.33e-2];
cell2Iter0Perf     = [8.12e-2, 6.85e-2, 5.91e-2, 5.80e-2, 6.10e-2, 6.77e-2, 8.16e-2, 9.78e-2];
cell2Iter1Perf     = [5.71e-2, 3.73e-2, 2.30e-2, 1.38e-2, 7.58e-3, 4.38e-3, 2.44e-3, 1.41e-3];
cell2Iter2Perf     = [5.67e-2, 3.82e-2, 2.28e-2, 1.31e-2, 7.03e-3, 3.55e-3, 1.71e-3, 8.49e-4];
cell2Iter2HardPerf = [5.51e-2, 3.61e-2, 2.05e-2, 1.01e-2, 4.60e-3, 1.80e-3, 6.71e-4, 2.46e-4];
thllr = [2.5, 2.5, 2.4, 2.3, 2.1, 1.8, 1.4, 0.9];
semilogy(dB, mmselimit, '-^', dB, zflimit, '-s', ...
        dB, cell2Iter0Perf, '-d', dB, cell2Iter1Perf, '-*', dB, cell2Iter2Perf, '-h', ...
        dB, cell2Iter2HardPerf, '-p', 'LineWidth', 1.5);
legend('无干扰场景, MMSE', '无干扰场景, ZF', ...
       '不进行干扰消除', '常规迭代1次', '常规迭代2次', '常规迭代两次，部分硬判1次');
set(gca, 'XTick', -10 : 2 : 4);
xlabel('SNR (dB)');
ylabel('BER');
axis([-10, 4, 1e-4, 1]);
grid on;

%% converge
% L = 2, M = K = 100, perfect channel estimation
asimu = [7.12, 7.12; 2.05, 2.05; 1.74, 1.74; 1.70, 1.71; 1.70, 1.70] - 1;
acalc = [7.13, 7.12; 2.04, 2.05; 1.70, 1.72; 1.70, 1.70; 1.70, 1.70] - 1;
asdb = 10 * log10(1 ./ asimu);
acdb = 10 * log10(1 ./ acalc);
xs  = [asdb(1, 1), asdb(2, 1), asdb(2, 1), asdb(4, 1), asdb(4, 1)];
ys  = [asdb(1, 2), asdb(1, 2), asdb(3, 2), asdb(3, 2), asdb(5, 2)];
xsp = [asdb(1, 1), asdb(1, 1), asdb(3, 1), asdb(3, 1), asdb(5, 1)];
ysp = [asdb(1, 2), asdb(2, 2), asdb(2, 2), asdb(4, 2), asdb(4, 2)];

xc  = [acdb(1, 1), acdb(2, 1), acdb(2, 1), acdb(4, 1), acdb(4, 1)];
yc  = [acdb(1, 2), acdb(1, 2), acdb(3, 2), acdb(3, 2), acdb(5, 2)];
xcp = [acdb(1, 1), acdb(1, 1), acdb(3, 1), acdb(3, 1), acdb(5, 1)];
ycp = [acdb(1, 2), acdb(2, 2), acdb(2, 2), acdb(4, 2), acdb(4, 2)];

plot(xs, ys, '-h', xsp, ysp, '-o', xc, yc, '-^', xcp, ycp, '-s', 'LineWidth', 1.5);
grid on;
legend('基站1，理论平均值', '基站2，理论平均值', '基站1，仿真数据', '基站2，仿真数据');
xlabel('SIR (dB)');
ylabel('SIR (dB)');


% L = 7, M = 100
% users = 76 : 4 : 100; % perfect channel estimation
users = 76 : 4 : 100; % imperfect channel estimation
Iter0M4 = [3.09e-2, 3.90e-2, 4.52e-2, 5.43e-2, 6.33e-2, 7.06e-2, 7.84e-2];
Iter1M4 = [4.79e-3, 7.42e-3, 1.19e-2, 1.75e-2, 2.40e-2, 3.13e-2, 3.78e-2];
Iter2M4 = [2.18e-3, 4.29e-3, 7.63e-3, 1.21e-2, 1.77e-2, 2.50e-2, 3.20e-2];
Iter0M0 = [2.23e-2, 3.00e-2, 3.78e-2, 4.58e-2, 5.36e-2, 6.16e-2, 7.06e-2];
Iter0M1 = [7.76e-4, 1.97e-3, 4.05e-3, 8.06e-3, 1.29e-2, 1.90e-2, 2.66e-2];
Iter0M2 = [1.13e-4, 3.95e-4, 1.26e-3, 3.34e-3, 7.08e-3, 1.26e-2, 1.96e-2];
semilogy(users, Iter0M4, '-^', users, Iter1M4, '-s', users, Iter2M4, '-d', users, Iter0M0, '-o', users, Iter0M1, '-v', users, Iter0M2, '-*', 'LineWidth', 1.5);
legend('SNR=0dB, #0', 'SNR=0dB, #1', 'SNR=0dB, #2', 'SNR=4dB, #0', 'SNR=4dB, #1', 'SNR=4dB, #2');
set(gca, 'XTick', 76 : 4 : 100);
%set(gca,'XTickLabel', {'76', '80', '84', '88', '92', '96', '100'});
xlabel('每小区用户数目(K)');
ylabel('BER');
axis([76, 100, 1e-4, 1e-1]);
grid on;


% Warning, obsolete
% SNR = 0dB,
users = 20 : 10 : 100;
free   = [-13.22, -14.91, -16.13, -17.08, -17.85, -18.51, -19.08, -19.59, -20.04];
cell2  = [-10.50, -11.24, -11.64, -11.91, -12.12, -12.31, -12.39, -12.43, -12.48];
cell3  = [ -8.99,  -9.53,  -9.77,  -9.93, -10.14, -10.20, -10.25, -10.33, -10.35];
%ocell2 = [-10.50, -11.15, -11.56, -11.88, -11.99, -12.23, -12.32, -12.40, -12.51] + 0.3;
%ocell3 = [ -8.94,  -9.32,  -9.70,  -9.99, -10.18, -10.25, -10.27, -10.35, -10.36] + 0.5;
crocell2est = [-15.14, -15.29, -15.48, -15.55, -15.56, -15.64, -15.66, -15.69, -15.76];
crocell = ones(1, 9) * -11.96;
%plot(users, free, '-^', users, cell2, '-d', users, cell3, '-s', users, crocell2est, '-v', users, crocell, 'LineWidth', 1.5);
plot(users, cell2, '-d', users, cell3, '-s', users, crocell2est, '-v', users, crocell, 'LineWidth', 1.5);
legend('Interference free scenario', 'L=2', 'L=3');
xlabel('Number of Users (K)');
ylabel('Mean Estimation Error (dB)');
grid on;


% K = 30, M = 40
dB = [-8 : 2 : 2];
oneCellImpR1= [1.01e-1, 6.30e-2, 3.48e-2, 1.53e-2, 5.42e-3, 1.21e-3];
oneCellImpR2= [7.89e-2, 4.59e-2, 2.21e-2, 8.30e-3, 2.24e-3, 4.07e-4];
oneCellImpR3= [6.91e-2, 3.81e-2, 1.78e-2, 6.78e-3, 1.56e-3, 2.61e-4];
oneCellImp  = [7.81e-2, 4.50e-2, 2.19e-2, 8.39e-3, 2.44e-3, 3.55e-4];
twoCellImp0 = [7.98e-2, 5.52e-2, 3.71e-2, 2.69e-2, 2.03e-2, 1.64e-2]; % L = 2, F = 1/2, # = 0
twoCellImp1 = [6.87e-2, 3.99e-2, 1.94e-2, 8.43e-3, 2.46e-3, 4.88e-4]; % L = 2, F = 1/2, # = 1
thrCellImp0 = [9.54e-2, 7.28e-2, 5.86e-2, 5.12e-2, 4.53e-2, 4.49e-2]; % L = 3, F = 1/3, # = 0
thrCellImp1 = [7.36e-2, 4.43e-2, 2.40e-2, 1.05e-2, 4.42e-3, 1.21e-3]; % L = 3, F = 1/3, # = 1
thrCellImp2 = [7.39e-2, 4.57e-2, 2.37e-2, 1.07e-2, 3.66e-3, 9.78e-4]; % L = 3, F = 1/3, # = 2
forCellImp0 = [9.75e-2, 8.41e-2, 6.89e-2, 6.30e-2, 6.29e-2, 6.29e-2]; % L = 4, F = 1/4, # = 0
forCellImp2 = [7.37e-2, 4.77e-2, 2.53e-2, 1.19e-2, 4.40e-3, 1.32e-3]; % L = 4, F = 1/4, # = 2
semilogy(dB, oneCellImpR1, '-d', dB, oneCellImpR2, '-h', ...
         dB, twoCellImp0,  '-s', dB, twoCellImp1,  '->', ...
         dB, thrCellImp0,  '-o', dB, thrCellImp2,  '-v', ...
         dB, forCellImp0,  '-p', dB, forCellImp2,  '-^', 'LineWidth', 1.5);
legend('L = 1, F = 1', 'L = 1, F = 1/2', ...
       'L = 2, F = 1/2, # = 0', 'L = 2, F = 1/2, # = 1', ...
       'L = 3, F = 1/3, # = 0', 'L = 3, F = 1/3, # = 2', ...
       'L = 4, F = 1/4, # = 0', 'L = 4, F = 1/4, # = 2');
set(gca, 'XTick', -8 : 0.2 : 2);
xlabel('SNR (dB)');
ylabel('BER');
axis([-8 2 1e-4 1.1e-1]);
grid on;


% M = 40, K = 30
dB = [-8 : 2 : 2];
oneCellImpM40K30 = [7.81e-2, 4.50e-2, 2.19e-2, 8.39e-3, 2.44e-3, 3.55e-4];
fourCellImp0R3   = [1.07e-1, 8.88e-2, 7.32e-2, 6.30e-2, 6.12e-2, 5.95e-2];
fourCellImp2R3   = [7.92e-2, 5.08e-2, 2.80e-2, 1.37e-2, 5.53e-3, 1.64e-3];
% M = 20, K = 10
oneCellImpM20R1  = [1.45e-1, 9.39e-2, 5.64e-2, 2.75e-2, 1.00e-2, 2.56e-3];
oneCellImpM20K10 = [9.12e-2, 5.17e-2, 2.46e-2, 8.31e-3, 1.91e-3, 2.49e-4];
oneCellPerM20K10 = [7.36e-2, 3.94e-2, 1.76e-2, 5.16e-3, 1.06e-3, 1.28e-4];

s7CellImp0R7     = [1.17e-1, 8.37e-2, 6.00e-2, 4.60e-2, 3.65e-2, 2.85e-2];
s7CellImp2R7     = [1.04e-1, 6.42e-2, 3.30e-2, 1.50e-2, 5.19e-3, 1.24e-3];
s7CellPerImp2    = [8.27e-2, 4.69e-2, 2.22e-2, 7.73e-3, 1.91e-3, 3.01e-4];

n19CellImp0R7    = [1.21e-1, 9.13e-2, 7.31e-2, 6.31e-2, 5.27e-2, 4.46e-2];
n19CellImp2R7    = [1.09e-1, 6.48e-2, 3.41e-2, 1.63e-2, 5.78e-3, 1.51e-3];
semilogy(dB, oneCellImpM20R1, '-d', dB, oneCellImpM20K10, '-h', dB, oneCellPerM20K10, '-s', ...
         dB, s7CellImp0R7, '->', dB, s7CellImp2R7, '-o', dB, s7CellPerImp2, '-v', ...
         dB, n19CellImp0R7, '-p', dB, n19CellImp2R7, '-^', 'LineWidth', 1.5);
legend('L = 1,  F = 1', 'L = 1,  F = 1/7', 'L = 1, 完美信道估计', ...
       'L = 7,  F = 1/7, # = 0', 'L = 7,  F = 1/7, # = 2', 'L = 7, 完美信道估计', ...
       'L = 19, F = 1/7, # = 0', 'L = 19, F = 1/7, # = 2');
set(gca, 'XTick', -8 : 0.2 : 2);
xlabel('SNR (dB)');
ylabel('BER');
axis([-8 2 1e-4 1.5e-1]);
grid on;



