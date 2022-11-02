function [Y, Thetas] = transmit_pilot(RIS_conf, BS_conf, f, G, Np)
    N = RIS_conf.N;  M = BS_conf.M;
    Ep = BS_conf.Ep; 
    pilot_power = Ep/Np*BS_conf.Pt_UE;
    Thetas = generate_RIS_codebook(N, Np, 'random');
    sigma_noise = BS_conf.sigma_noise;
    
    % Channel estimation scheme: Just estimate fk and G together, for each user.
    HT = (G.') * diag(conj(f));
    tmp = sqrt(pilot_power)*HT*Thetas;
    Y = tmp + sigma_noise * (randn([M, Np])+1j*randn([M, Np])/sqrt(2));
end

