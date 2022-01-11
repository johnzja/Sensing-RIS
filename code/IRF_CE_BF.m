function [P_recv, Rate] = IRF_CE_BF(N, M, RIS_conf, BS_conf, f, G, Np)
    % Fix the BS precoding to steer the beam at RIS.
    BS_pos_sph = BS_conf.BS_pos_sph;
    theta_BS = BS_pos_sph(2); psi_BS = BS_pos_sph(3);
    My = M(1); Mz = M(2);
    
    w = kron(exp(-1j*pi*sin(psi_BS)*(0:Mz-1)), exp(-1j*pi*cos(psi_BS)*sin(theta_BS)*(0:My-1)));
    
    % Generate & Measure the interferential random field (IRF) on each RIS
    % element.
    
    




end