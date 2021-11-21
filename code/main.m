%%% Sensing RIS
addpath('code/');

% DEMO version.
% Setup simulation parameters.
Ts = 1e-3;
alpha = 1;
beta = 0.5;

L = 2^8;
f_psi = 1/Ts;

A = 1;

sigma_v = 0.4;
sigma_zeta = 0.05;

randn('seed', 0);
rand('seed', 0);

N_exp = 2000;
MSE_arr_LS = zeros(N_exp, 1);
MSE_arr_Newton = zeros(N_exp, 1);

for idx = 1:N_exp
% Generate power-sensor signal.
P = zeros(L, 1);
varphi = 2*pi*rand();   % uniform (0, 2pi).

for l = 1:L
    v = (sigma_v*randn() + 1j*sigma_v*randn())/sqrt(2);
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
MSE_arr_LS(idx) = (delta - round(delta))^2;

% Study the correctness of the derivative.
% logL = zeros(L, 1);
% dlogL = zeros(L, 1);
% d2logL = zeros(L, 1);
% phase_arr = 2*pi*(0:L-1).'/L;
% for l = 1:L
%     vphi_hat = phase_arr(l);
%     [a, b, c] = calc_likelihood(P, alpha, beta, A, 2*pi*f_psi*(0:L-1).'*Ts/L, vphi_hat, sigma_v);
%     logL(l) = a;
%     dlogL(l) = b;
%     d2logL(l) = c;
% end
% figure(1);
% plot(phase_arr/(pi), logL);
% xlabel('x\pi radians');
% 
% figure(2);
% plot(phase_arr/pi, dlogL);
% xlabel('x\pi radians');
% 
% figure(3);
% plot(phase_arr/pi, d2logL);
% xlabel('x\pi radians');

% Estimate varphi by Newton-Raphson method.
for k = 1:3
    [logL, dlogL, d2logL] = calc_likelihood(P, alpha, beta, A, 2*pi*f_psi*(0:L-1).'*Ts/L, varphi_hat, sigma_v);
    varphi_hat = varphi_hat - dlogL/d2logL;
end

delta = (varphi_hat - varphi)/(2*pi);
MSE_arr_Newton(idx) = (delta - round(delta))^2;

end

fprintf('Sim Complete.\n');
fprintf('std var LS  \t = %f degs\n', sqrt(mean(MSE_arr_LS))*360);
fprintf('std var Newton\t = %f degs\n', sqrt(mean(MSE_arr_Newton))*360);


