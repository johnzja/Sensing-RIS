function CRLB = get_CRLB(SensingRIS_param, varphi)
% Calculate the CRLB (in radians). 
% Parameter extraction.
    A = SensingRIS_param.A;
    alpha = SensingRIS_param.alpha;
    beta = SensingRIS_param.beta;
    sigma_v = SensingRIS_param.sigma_v;
    psi_arr = SensingRIS_param.psi_arr;

    sine_squared = sin(psi_arr + varphi).^2;
    a = A*sigma_v^2;
    lambda = A*(alpha^2 + beta^2 + 2*alpha*beta*cos(psi_arr + varphi));
    gamma = lambda / a;     % Interference-SNR.
    g = sqrt(pi./gamma) .* exp(-gamma/2) .* ((1+1./gamma).*besseli(0, gamma/2) + besseli(1, gamma/2))/4;
    
    tmp = a./lambda - g;
    tmp = tmp .* (tmp>0);
    CRLB = 1/(sine_squared.' * tmp) * (sigma_v^4)/(4*alpha^2*beta^2);
end

