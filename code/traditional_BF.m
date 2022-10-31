function [P_recv, Rate] = traditional_BF(RIS_conf, BS_conf, f, G)
% P - transmit power
% f - N*1
% G - N*M
% w - M*1
% theta - N*1
    [N, ~] = size(G);
    assert(N == RIS_conf.N);
    % Perform beamforming based on H_hat. The model is: y_UE = theta.' * (H).' *w*s + n
    theta = exp(1j*2*pi*rand(N, 1));
    iter = 10;
    Pt_BS = BS_conf.Pt_BS;
    P_noise = BS_conf.sigma_noise^2;
    
    % Assume that s=1.
    objective = zeros(2*iter, 1);
    threshold = 0.01;
    H = G.'*diag(conj(f));
    
    for idx = 1:iter
        % update w
        w = (theta' * H').';
        w = sqrt(Pt_BS)*w/norm(w);
        objective(2*idx-1) = abs(theta.'*H.'*w)^2;
        % update theta
        theta = exp(-1j*angle(H.'*w));
        t = theta.'*H.'*w;
        objective(2*idx) = abs(t)^2;
        if abs(objective(2*idx)-objective(2*idx-1))/objective(2*idx) < threshold
            P_recv = abs(theta.'*H.'*w)^2;         % Incorporate the true channel.
            Rate = log2(1+P_recv/P_noise);
            return;
        end
    end
    P_recv = abs(theta.'*H.'*w)^2;
    Rate = log2(1+P_recv/P_noise);  
end

