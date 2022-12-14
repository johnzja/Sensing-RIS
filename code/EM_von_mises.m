function varphi = EM_von_mises(P, SensingRIS_param, N_iters, FFT_first)
% Calculate the estimated varphi using von-mises EM algorithm. 
% Parameter extraction.
    A = SensingRIS_param.A;
    alpha = SensingRIS_param.alpha;
    beta = SensingRIS_param.beta;
    sigma_v = SensingRIS_param.sigma_v;
    psi_arr = SensingRIS_param.psi_arr;

    if ~exist("N_iters", "var")
        N_iters = 4;
    end
    
    if ~exist("FFT_first", "var")
        FFT_first = true;
    end

    % Step 1: Obtain coarse estimation for varphi.

    if FFT_first
        p = fft(P);
        varphi_est = angle(p(2));
        varphi_kappa = 1;       % Adjustable parameter.
    else
        varphi_est = 0; 
        varphi_kappa = 0;       % Without any prior knowledge. 
    end
    
    s = sqrt(P/A);  % estimator for abs(alpha + beta*exp(j??)*exp(j??(t))).
    L = length(P);
    assert(L == length(psi_arr));
    
    thetas_mu = zeros(L, 1);
    thetas_kappa = zeros(L, 1);
    
    % Step 2: Update estimation by von-mises Bayesian algorithm.
    for iter = 1:N_iters
        % Estimate each theta[l] by Bayesian rule.
        mus = alpha + beta*exp(1j*varphi_est)*exp(1j*psi_arr);
        for idx = 1:L
            thetas_kappa(idx) = s(idx) / (sigma_v^2/2) * abs(mus(idx));
            thetas_mu(idx) = angle(mus(idx));
        end

        % Estimate varphi by Bayesian rule.
        obs_points = s .* exp(1j*thetas_mu) - alpha;    % subtract the mean.
        % scatter(real(obs_points), imag(obs_points));
        cv_varphi = varphi_kappa*exp(1j*varphi_est) + beta*sum(obs_points.*exp(-1j*psi_arr))/(sigma_v^2/2);

        varphi_est = angle(cv_varphi);
        % varphi_kappa = abs(cv_varphi);
    end
    varphi = varphi_est;
end