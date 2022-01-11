function [P_recv, Rate] = traditional_CE_BF(RIS_conf, BS_conf, f, G, Np)
    % Extract necessary parameters.
    N = RIS_conf.N;
    F = dftmtx(N);
    Thetas = F(:, 1:Np);    
    sigma_noise = BS_conf.sigma_noise;
    P_noise = BS_conf.sigma_noise^2;
    M = BS_conf.M;
    
    % Channel estimation scheme: Just estimate fk and G together, for each user.
    H = (G.') * diag(conj(f));
    tmp = sqrt(BS_conf.Pt_UE)*H*Thetas;
    Y = tmp + sigma_noise * (randn([M, Np])+1j*randn([M, Np])/sqrt(2));
    H_hat = Y*Thetas'/Np/sqrt(BS_conf.Pt_UE);
    
    % Perform beamforming based on H_hat. The model is: y_UE = theta.' * (H).' *w*s + n
    theta = exp(1j*2*pi*rand(N, 1));
    iter = 10;
    Pt_BS = BS_conf.Pt_BS;
    
    % Assume that s=1.
    objective = zeros(2*iter, 1);
    threshold = 0.01;
    for idx = 1:iter
        % update w
        w = (theta' * H_hat').';
        w = sqrt(Pt_BS)*w/norm(w);
        objective(2*idx-1) = abs(theta.'*H_hat.'*w)^2;
        % update theta
        theta = exp(-1j*angle(H_hat.'*w));
        t = theta.'*H_hat.'*w;
        objective(2*idx) = abs(t)^2;
        if abs(objective(2*idx)-objective(2*idx-1))/objective(2*idx) < threshold
            P_recv = abs(theta.'*H.'*w)^2;         % Incorporate the true channel.
            Rate = log2(1+P_recv/P_noise);
            return;
        end
    end
    P_recv = abs(theta.'*H.'*w)^2;
    Rate = log2(1+P_recv/P_noise);  % In fact, this is the spectral efficiency.
end