function [w, theta] = RIS_precode(RIS_conf, HT, initial_theta)
    N = RIS_conf.N; 
    
    % Perform beamforming based on HT_hat. The model is: y_UE = theta.' * (HT).' *w*s + n
    if ~exist('initial_theta', 'var')
        theta = exp(1j*2*pi*rand(N, 1));
    else
        theta = initial_theta;
    end
    iter = 10;
    
    % Assume that s=1.
    objective = zeros(2*iter, 1);
    threshold = 0.01;
    % HT = G.'*diag(conj(f_hat));
    
    for idx = 1:iter
        % update w
        w = (theta' * HT').';
        w = w/norm(w);              % Assume unit-length. 
        objective(2*idx-1) = abs(theta.'*HT.'*w)^2;
        % update theta
        theta = exp(-1j*angle(HT.'*w));
        t = theta.'*HT.'*w;
        objective(2*idx) = abs(t)^2;
        if abs(objective(2*idx)-objective(2*idx-1))/objective(2*idx) < threshold
            return;
        end
    end
    
end

