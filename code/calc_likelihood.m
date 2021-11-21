function [L,dL,d2L] = calc_likelihood(P, alpha, beta, A, psi_arr, varphi, sigma_v)
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

