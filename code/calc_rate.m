function R = calc_rate(BS_conf, G, f, w, theta, Np)
% Ensure that |w|=1. 
    HT = G.'*diag(conj(f)); % true cascaded channel.
    gamma = BS_conf.Pt_BS*abs(theta.'*(HT.')*w)^2/BS_conf.sigma_noise^2;  % Receive SNR.
    BL = BS_conf.BL;
    R = (BL-Np)/BL*log2(1+gamma);
end
