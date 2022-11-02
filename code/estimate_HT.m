function [HT_hat, nmse] = estimate_HT(RIS_conf, BS_conf, Y, Thetas, f, G, Np, alg_est)
    N = RIS_conf.N;
    Ep = BS_conf.Ep; 
    lambda = RIS_conf.lambda;
    
    pilot_power = Ep/Np*BS_conf.Pt_UE;
    P_noise = BS_conf.sigma_noise^2;
    
    % Channel estimation scheme: Just estimate fk and G together, for each user.
    HT = (G.') * diag(conj(f));
    
    switch alg_est
        case 'MF'
            HT_hat = Y*Thetas'/Np/sqrt(pilot_power);    % Matched filter. 
        case 'LS'
            FFH = Thetas*Thetas';
            HT_hat = Y*Thetas'/(FFH)/sqrt(pilot_power);
        case 'MMSE'
            FFH = Thetas*Thetas';
            sigma2_h = (lambda/(4*pi*100))^4;
            HT_hat = sqrt(pilot_power)*Y*Thetas'/(FFH*pilot_power + P_noise/sigma2_h*eye(N));
        otherwise
            error('parameter alg_est error\n');
    end

    % Calculate estimation noise here, in NMSE. 
    nmse = pow2db((norm(HT_hat-HT, 'fro')/norm(HT,'fro'))^2);
    
%     % The model is: y_UE = theta.' * (H).' *w*s + n
%     theta = exp(1j*2*pi*rand(N, 1));
%     iter = 10;
%     
%     % Assume that E|s|^2 = 1.
%     objective = zeros(2*iter, 1);
%     threshold = 0.01;
%     
%     for idx = 1:iter
%         % update w
%         w = (theta' * HT_hat').';
%         w = sqrt(data_power)*w/norm(w);             % Transmit power integrated. 
%         objective(2*idx-1) = abs(theta.'*HT_hat.'*w)^2;
%         % update theta
%         theta = exp(-1j*angle(HT_hat.'*w));
%         t = theta.'*HT_hat.'*w;
%         objective(2*idx) = abs(t)^2;
%         if abs(objective(2*idx)-objective(2*idx-1))/objective(2*idx) < threshold
%             P_recv = abs(theta.'*HT.'*w)^2;         % Incorporate the true channel.
%             Rate = (BL-Np)/BL*log2(1+P_recv/P_noise);
%             return;
%         end
%     end
%     
%     P_recv = abs(theta.'*HT.'*w)^2;
%     Rate = (BL-Np)/BL*log2(1+P_recv/P_noise);  % In fact, this is the spectral efficiency.

end

