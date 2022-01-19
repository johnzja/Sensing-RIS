function CRLB = get_precise_CRLB(SensingRIS_param, varphi)
% Calculate the CRLB (in radians). 
% Parameter extraction.
% Load the data.
persistent g_values;
persistent gamma_range;
if ~exist('g_values', 'var') || isempty(g_values)
    load data/preciseCRLB_data.mat;
    gamma_range = 0.1:0.01:5;
    g_values = Expression1;
end

    A = SensingRIS_param.A;
    alpha = SensingRIS_param.alpha;
    beta = SensingRIS_param.beta;
    sigma_v = SensingRIS_param.sigma_v;
    psi_arr = SensingRIS_param.psi_arr;

    sine_squared = sin(psi_arr + varphi).^2;
    a = A*sigma_v^2;
    lambda = A*(alpha^2 + beta^2 + 2*alpha*beta*cos(psi_arr + varphi));
    gamma = lambda / a;     % Interference-SNR.
    % Replace the g function with the exact value.
    tmp = zeros(size(gamma));
    for idx = 1:length(gamma)
        if 0<=gamma(idx) && gamma(idx) <= gamma_range(end)
            tmp(idx) = interp1(gamma_range, g_values, gamma(idx), 'spline', 'extrap');
        else
            tmp(idx) = 1/gamma(idx) - sqrt(pi/gamma(idx)) * exp(-gamma(idx)/2) * ((1+1/gamma(idx))*besseli(0,gamma(idx)/2) + besseli(1,gamma(idx)/2))/4;
        end
    end
    % g = sqrt(pi./gamma) .* exp(-gamma/2) .* ((1+1./gamma).*besseli(0, gamma/2) + besseli(1, gamma/2))/4;
    CRLB = 1/(sine_squared.' * tmp) * (sigma_v^4)/(4*alpha^2*beta^2);
end
