%%% Sensing RIS
% Setup simulation parameters.
rng(0);

Ts = 1e-3;
alpha = 1;
A = 1;
L = 2^8;
f_psi = 1/Ts;       % Frequency offset w.r.t carrier-freq.
sigma_zeta  = 0.05;
N_exp = 4000;       % Number of numerical experiments.

% Scan parameters.
K = 0.6;
gamma_bar = 2;
N_scan = 9;        % # Scan points.

% K_arr = linspace(0.2,0.5,N_scan);
gamma_bar_arr = linspace(1, 5, N_scan);
MSE_arr_LS = zeros(N_exp, N_scan);
MSE_arr_VM = zeros(N_exp, N_scan);
MSE_arr_Newton = zeros(N_exp, N_scan);
CRLB = zeros(1, N_scan);
CRLB_precise = zeros(1, N_scan);
    
for idx_scan = 1:N_scan
    % K = K_arr(idx_scan);
    gamma_bar = gamma_bar_arr(idx_scan);
    beta    = (1-sqrt(1-K^2))/K;
    sigma_v = sqrt((alpha^2+beta^2)/gamma_bar);

    SensingRIS_param.alpha      = alpha;
    SensingRIS_param.beta       = beta;
    SensingRIS_param.A          = A;
    SensingRIS_param.L          = L;
    SensingRIS_param.Ts         = Ts;
    SensingRIS_param.f_psi      = f_psi;
    SensingRIS_param.psi_arr    = 2*pi*f_psi*(0:L-1).'*Ts/L;
    SensingRIS_param.sigma_v    = sigma_v;

    % Generate power-sensor signal.
    P = zeros(L, 1);
    varphi = 2*pi*rand();   % uniform (0, 2pi).
    CRLB(idx_scan) = get_CRLB(SensingRIS_param, varphi);
    CRLB_precise(idx_scan) = get_precise_CRLB(SensingRIS_param, varphi);
    
    for idx = 1:N_exp
        % Generate power signals.
        for l = 1:L
            v = (randn() + 1j*randn())*sigma_v/sqrt(2);
            P(l) = A*abs(alpha + beta*exp(1j*(2*pi*f_psi*(l-1)*Ts/L + varphi)) + v)^2;
            if sigma_zeta > 0
                P(l) = P(l) + sigma_zeta * randn();
                if P(l)<0
                    P(l) = 1e-5;
                end
            end
        end

        % LS method with FFT.
        p = fft(P);
        varphi_hat = angle(p(2));

        delta = (varphi_hat - varphi)/(2*pi);
        MSE_arr_LS(idx, idx_scan) = ((delta - round(delta))*(2*pi))^2;

        % von-Mises method
        varphi_hat = EM_von_mises(P, SensingRIS_param);
        delta = (varphi_hat - varphi)/(2*pi);
        MSE_arr_VM(idx, idx_scan) = ((delta - round(delta))*(2*pi))^2;

        % Estimate varphi by Newton-Raphson method.
        for k = 1:4
            [logL, dlogL, d2logL] = calc_likelihood(P, varphi_hat, SensingRIS_param);
            varphi_hat = varphi_hat - dlogL/d2logL;
        end

        delta = (varphi_hat - varphi)/(2*pi);
        MSE_arr_Newton(idx, idx_scan) = ((delta - round(delta))*(2*pi))^2;
    end

    fprintf('Sim Complete for gamma_bar=%f, K=%f.\n', gamma_bar, K);
    fprintf('std var LS  \t = %f rads\n', sqrt(mean(MSE_arr_LS(:,idx_scan))));
    fprintf('std var VM  \t = %f rads\n', sqrt(mean(MSE_arr_VM(:,idx_scan))));
    fprintf('std var Newton\t = %f rads\n', sqrt(mean(MSE_arr_Newton(:,idx_scan))));
    fprintf('std var CRLB\t = %f rads\n', sqrt(CRLB_precise(idx_scan)));
    fprintf('----------------------------------------------\n');
end
%% Plot the results.

MSE_LS = mean(MSE_arr_LS);
MSE_VM = mean(MSE_arr_VM);
MSE_Newton = mean(MSE_arr_Newton);

set(0,'DefaultLineMarkerSize',4);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth',1.4);
set(0,'defaultfigurecolor','w');
figure('color',[1 1 1]); hold on;
plot(gamma_bar_arr, MSE_LS, 'bp-');
plot(gamma_bar_arr, MSE_VM, 'gs-');
plot(gamma_bar_arr, MSE_Newton, 'ro-');
plot(gamma_bar_arr, CRLB_precise, 'ko-.');
plot(gamma_bar_arr, CRLB, 'color', [228,0,127]/255, 'LineStyle', '--', 'marker', 'o');

set(gca,'FontName','Times New Roman');
grid on; box on;
legend('LS', 'VM-EM', 'Newton-ML', 'CRLB', 'CRLB-approx');
xlabel('$\bar{\gamma}$', 'interpreter', 'latex');
% xlabel('$K$', 'interpreter', 'latex');
ylabel('E$\left(\hat{\varphi}-\varphi\right)^2$', 'interpreter', 'latex');

