function [L,dL,d2L] = calc_likelihood(P, varphi, SensingRIS_param)
% Calculate the log-likelihood function of observed power sequence P under phase assumption varphi.
% Parameter extraction.
    A = SensingRIS_param.A;
    alpha = SensingRIS_param.alpha;
    beta = SensingRIS_param.beta;
    sigma_v = SensingRIS_param.sigma_v;
    psi_arr = SensingRIS_param.psi_arr;
    
    lambda = A*(alpha^2 + beta^2 + 2*alpha*beta*cos(psi_arr + varphi));
    z = sqrt(P .* lambda)/(A*sigma_v^2/2);
    len = length(P);
    assert(len == length(psi_arr));
    
    L = -len * log(A*sigma_v^2) - sum(P+lambda)/(A*sigma_v^2) + sum(log(besseli(0,z)));
    
    tmp = (2*alpha*beta)/(sigma_v^2);
    R = besseli(1,z) ./ besseli(0,z);
    P_over_lambda = P./lambda;
    sqrt_P_over_lambda = sqrt(P_over_lambda);
    
    dL = tmp * (sin(psi_arr + varphi).')*(1-(R.*sqrt_P_over_lambda));
    
    d2L = tmp * (cos(psi_arr+varphi).') * (1-(R.*sqrt_P_over_lambda));
    d2L = d2L + tmp^2 * (sin(psi_arr+varphi).^2).'*(P_over_lambda.*(1-R.^2-2*R./z));
end

