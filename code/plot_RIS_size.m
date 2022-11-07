
N_arr = [10*5, 10*10, 20*10, 30*10, 40*10, 40*20];
Rate_IRF = zeros(6, 1);
Rate_Oracle = zeros(6, 1);
Rate_MMSE_f = zeros(6, 1);
Rate_LMMSE_H = zeros(6, 1);

%% 10x5
load('data/IRF_sim-10x5.mat', 'SE');

Rate_IRF(1) = SE(5, 6);
Rate_Oracle(1) = SE(5, 7);
Rate_MMSE_f(1) = SE(5, 5);
Rate_LMMSE_H(1) = SE(5, 4); 

%% 10x10
load('data/IRF_sim-10x10.mat', 'SE');

Rate_IRF(2) = SE(5, 6);
Rate_Oracle(2) = SE(5, 7);
Rate_MMSE_f(2) = SE(5, 5);
Rate_LMMSE_H(2) = SE(5, 4); 

%% 20x10
load('data/IRF_sim-20x10.mat', 'SE');

Rate_IRF(3) = SE(5, 6);
Rate_Oracle(3) = SE(5, 7);
Rate_MMSE_f(3) = SE(5, 5);
Rate_LMMSE_H(3) = SE(5, 4); 

%% 30x10
load('data/IRF_sim-30x10.mat', 'SE');

Rate_IRF(4) = SE(5, 6);
Rate_Oracle(4) = SE(5, 7);
Rate_MMSE_f(4) = SE(5, 5);
Rate_LMMSE_H(4) = SE(5, 4); 

%% 40x10
load('data/IRF_sim-40x10.mat', 'SE');

Rate_IRF(5) = SE(5, 6);
Rate_Oracle(5) = SE(5, 7);
Rate_MMSE_f(5) = SE(5, 5);
Rate_LMMSE_H(5) = SE(5, 4); 

%% 40x20
load('data/IRF_sim-40x20.mat', 'SE');

Rate_IRF(6) = SE(5, 6);
Rate_Oracle(6) = SE(5, 7);
Rate_MMSE_f(6) = SE(5, 5);
Rate_LMMSE_H(6) = SE(5, 4); 

%% Plots.

set(0,'DefaultLineMarkerSize',6);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth',1.4);
set(0,'defaultfigurecolor','w');

figure('color',[1 1 1]); hold on;

plot(N_arr, Rate_Oracle, 'ko-.','MarkerSize',6, 'LineWidth', 2);                                % Oracle
plot(N_arr, Rate_IRF, 'ro-','MarkerSize',6);                                                    % IRF
plot(N_arr, Rate_MMSE_f, 'color', [0, 0, 1], 'LineStyle', '-', 'marker', 'x','MarkerSize',8);   % MMSE f
plot(N_arr, Rate_LMMSE_H, 'color', [1, 0, 1], 'LineStyle', '--', 'marker', 'p','MarkerSize',6); % LMMSE H

set(gca,'FontName','Times New Roman');
set(gca, 'xscale', 'log');
grid on; box on;
legend('Oracle', 'Proposed-IRF + vM-EM', 'MMSE f', 'LMMSE H');
xlabel('$N_{\rm RIS}$', 'interpreter', 'latex');
ylabel('$C_{\rm SI}$ (bps/Hz)', 'interpreter', 'latex');

