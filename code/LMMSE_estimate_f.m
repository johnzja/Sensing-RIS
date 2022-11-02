function [f_hat, Np] = LMMSE_estimate_f(RIS_conf, BS_conf, f, G, Ep, Ed)
    M = BS_conf.M; N = RIS_conf.N; BL = BS_conf.BL; 
    lambda = RIS_conf.lambda;
    Np = ceil(N/M);
    sigma_noise = BS_conf.sigma_noise;
    
    Thetas = generate_RIS_codebook(N, Np);
    
    pilot_power = Ep/Np*BS_conf.Pt_UE;      % uplink channel estimation. 
    data_power = Ed/(BL-Np)*BS_conf.Pt_BS;  % downlink data transfer. 
    
    
    GT = G.';
    A = zeros(M*Np, N);
    y_vec = zeros(M*Np, 1);
    f_star = conj(f);
    
    for idx = 1:Np
        TMP = GT*diag(Thetas(:, idx));
        n = (randn(M,1) + 1j*randn(M,1))/sqrt(2);
        A((idx-1)*M+1:idx*M, :) = TMP;
        y_vec((idx-1)*M+1:idx*M) = sqrt(pilot_power)*(TMP*f_star) + sigma_noise*n;
    end
    
    % Estimate f_star. 
    % LMMSE formula: x = H\theta+w => \hat{\theta} = (H'*H + 1/gamma*eye(p))\(H'*x);
    % gamma = sigma2_theta/sigma2_w = transmit_signal_power / receive_noise_power.
    gamma = (lambda/(4*pi*100))^2/sigma_noise^2;
    f_star_hat = (A'*A + eye(N)/gamma)\(A'*y_vec);
    f_hat = conj(f_star_hat);
end
