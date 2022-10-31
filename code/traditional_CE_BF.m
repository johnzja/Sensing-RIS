function [P_recv, Rate] = traditional_CE_BF(RIS_conf, BS_conf, f, G, Ep, Ed, Np)
    % Extract necessary parameters.
    N = RIS_conf.N;
    M = BS_conf.M;
    lambda = RIS_conf.lambda;
    BL = BS_conf.BL; 
    
    F = dftmtx(N);
    pilot_power = Ep/Np*BS_conf.Pt_UE;
    data_power = Ed/(BL-Np)*BS_conf.Pt_BS;  % downlink power. 
    
    if Np <= N
        Thetas = F(:, 1:Np);    
    else
        Thetas = zeros(N, Np);
        Thetas(:,1:N) = F;
        Thetas(:, N+1:end) = exp(1j*2*pi*rand([N, Np-N]));
    end
    
    sigma_noise = BS_conf.sigma_noise;
    P_noise = BS_conf.sigma_noise^2;
    
    % Channel estimation scheme: Just estimate fk and G together, for each user.
    HT = (G.') * diag(conj(f));
    tmp = sqrt(pilot_power)*HT*Thetas;
    Y = tmp + sigma_noise * (randn([M, Np])+1j*randn([M, Np])/sqrt(2));
    
    HT_hat_LS = Y*Thetas'/Np/sqrt(pilot_power);   % LS-CE (approx)  
    
    FFH = Thetas*Thetas';
    sigma2_h = (lambda/(4*pi*100))^4;
    HT_hat_LS_precise = Y*Thetas'/(FFH)/sqrt(pilot_power);
    HT_hat_MMSE = sqrt(pilot_power)*Y*Thetas'/(FFH*pilot_power + P_noise/sigma2_h*eye(N));

    % Calculate estimation noise here, in NMSE. 
    nmse_db_LS = pow2db((norm(HT_hat_LS-HT, 'fro')/norm(HT,'fro'))^2);
    nmse_db_LSp = pow2db((norm(HT_hat_LS_precise-HT, 'fro')/norm(HT,'fro'))^2);
    nmse_db_MMSE = pow2db((norm(HT_hat_MMSE-HT, 'fro')/norm(HT,'fro'))^2);
    
    % Perform beamforming based on the best HT_hat_MMSE. The model is: y_UE = theta.' * (H).' *w*s + n
    theta = exp(1j*2*pi*rand(N, 1));
    iter = 10;
    
    % Assume that E|s|^2 = 1.
    HT_hat = HT_hat_MMSE;
    objective = zeros(2*iter, 1);
    threshold = 0.01;
    
    for idx = 1:iter
        % update w
        w = (theta' * HT_hat').';
        w = sqrt(data_power)*w/norm(w);             % Transmit power integrated. 
        objective(2*idx-1) = abs(theta.'*HT_hat.'*w)^2;
        % update theta
        theta = exp(-1j*angle(HT_hat.'*w));
        t = theta.'*HT_hat.'*w;
        objective(2*idx) = abs(t)^2;
        if abs(objective(2*idx)-objective(2*idx-1))/objective(2*idx) < threshold
            P_recv = abs(theta.'*HT.'*w)^2;         % Incorporate the true channel.
            Rate = (BL-Np)/BL*log2(1+P_recv/P_noise);
            return;
        end
    end
    
    P_recv = abs(theta.'*HT.'*w)^2;
    Rate = (BL-Np)/BL*log2(1+P_recv/P_noise);  % In fact, this is the spectral efficiency.
end