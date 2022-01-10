function [P_recv, Rate] = traditional_channel_estimation_beamforming(N, M, RIS_conf, BS_conf, f, G, Np)
    Ny = N(1); Nz = N(2);
    My = M(1); Mz = M(2);
    N = Ny * Nz;
    
    F = dftmtx(N);
    Thetas = F(:, 1:Np);    
    
    % Channel estimation scheme: Just estimate fk and G together, for each user.
    PSD_noise = db2pow(-174-30);
    BW = 100e6;                 % System baseband BW = 100MHz.
    P_noise = PSD_noise * BW;   % Thermonoise for BS receiver equipped with M RF-chains.
    sigma_noise = sqrt(P_noise);
    
    H = (G.') * diag(conj(f));
    Y = H*Thetas + sigma_noise * randn([M, Np]);
    H_hat = Y*Thetas'/Np;
    
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
        t = theta.'*H_hat.'*w;
        theta = exp(-1j*angle(t));
        objective(2*idx) = abs(t)^2;
        if abs(objective(2*idx)-objective(2*idx-1))/objective(2*idx) < threshold
            P_recv = objective(2*idx);
            Rate = log2(1+P_recv/P_noise);
            return;
        end
    end
    P_recv = objective(2*iter);
    Rate = log2(1+P_recv/P_noise);
end