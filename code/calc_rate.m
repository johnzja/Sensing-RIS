function R = calc_rate(BS_conf, G, f, w, theta, Np)
% Ensure that |w|=1. 
    BL = BS_conf.BL;
    HT = G.'*diag(conj(f));                                                 % true cascaded channel.
    data_power = BS_conf.Ed/(BL-Np)*BS_conf.Pt_BS;                          % downlink data transfer. 
    gamma = data_power*abs(theta.'*(HT.')*w)^2/BS_conf.sigma_noise^2;       % Receive SNR.
    R = (BL-Np)/BL*log2(1+gamma);
end
