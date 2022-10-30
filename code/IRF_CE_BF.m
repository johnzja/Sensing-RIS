function [P_recv, Rate] = IRF_CE_BF(RIS_conf, BS_conf, f, G, Np, channel_type)
    % Fix the BS precoding to steer the beam at RIS.
    BS_pos_sph = BS_conf.BS_pos_sph;
    theta_BS = BS_pos_sph(2); psi_BS = BS_pos_sph(3);
    My = BS_conf.My; Mz = BS_conf.Mz;
    Ny = RIS_conf.Ny; Nz = RIS_conf.Nz; N_ris = RIS_conf.N;
    
    % Calculate the estimated optimal precoding directly from the BS-RIS geometries.
    if strcmp(channel_type, 'SV')
        w = kron(exp(1j*pi*sin(psi_BS)*(0:Mz-1)), exp(1j*pi*cos(psi_BS)*sin(theta_BS)*(0:My-1))).';
        w = sqrt(BS_conf.Pt_BS)*w/norm(w);
    elseif strcmp(channel_type, 'rayleigh')
        [~, ~, V] = svd(G);
        w = V(:,1);
        w = sqrt(BS_conf.Pt_BS)*w;
    end
    
    A = 1;
    L = 256;        % L samples within one OFDM symbol. 
	sigma_v = RIS_conf.sigma_v;
    
    psi_arr = 2*pi*(0:L-1).'/L;
    exp_psi_arr = exp(1j*psi_arr);
    SensingRIS_param.psi_arr    = psi_arr;
    SensingRIS_param.A          = A;
    SensingRIS_param.L          = L;
    SensingRIS_param.sigma_v    = sigma_v;
    
    % Generate & Measure the interferential random field (IRF) on each RIS
    % element.
    
    % Calculate the alpha-field.
    s_alpha = G*w;
    
    % Calculate the beta-field.
    s_beta = conj(f)*sqrt(BS_conf.Pt_UE);

    % Report the mean interferential SNR. 
    gamma_bar_arr = zeros(N_ris, 1);
    K_arr = zeros(N_ris, 1);

    % Obtain the interference random field (IRF) for each time slot.
    P_int = zeros(L, N_ris);
    for idx = 1:L
        P_int(idx,:) = abs(s_alpha.' + exp_psi_arr(idx)*s_beta.' + ...
            (randn([1,N_ris])+1j*randn([1,N_ris]))/sqrt(2)*sigma_v).^2;
    end
    
    % Estimate the phase difference varphi for each of the RIS elements.
    varphi = zeros(N_ris, 1);
    for n=1:N_ris
        SensingRIS_param.alpha	= abs(s_alpha(n)+sigma_v*(randn()+1j*randn())/sqrt(2*L));  
        % Noises should also be added here.
        SensingRIS_param.beta	= abs(s_beta(n)+sigma_v*(randn()+1j*randn())/sqrt(2*L));

        K_arr(n) = 2*(SensingRIS_param.alpha)*(SensingRIS_param.beta)/(SensingRIS_param.alpha^2+SensingRIS_param.beta^2);
        gamma_bar_arr(n) = (SensingRIS_param.alpha^2+SensingRIS_param.beta^2)/sigma_v^2;

        varphi(n) = EM_von_mises(P_int(:,n), SensingRIS_param);
    end
    
    K_mean = mean(K_arr);
    gamma_bar_mean = mean(gamma_bar_arr);
    % report the results. 
    % fprintf('K_mean = %f, \\gamma_bar_mean = %f\n', K_mean, gamma_bar_mean);
    
    % Calculate the (invariant) G.
    if strcmp(channel_type, 'SV')
        phase_G = kron(exp(1j*pi*sin(psi_BS)*(0:Nz-1)), exp(1j*pi*cos(psi_BS)*sin(theta_BS)*(0:Ny-1)));
        phase_BS = kron(exp(-1j*pi*sin(psi_BS)*(0:Mz-1)), exp(-1j*pi*cos(psi_BS)*sin(theta_BS)*(0:My-1)));
        G_phases_only = kron(phase_G.', phase_BS);  % size(G) = [N, M].

        % Calculate the best RIS phases theta_vec, directly on the RIS.
        arg_gnT_w = angle(G_phases_only*w);
    elseif strcmp(channel_type, 'rayleigh')
        arg_gnT_w = angle(G*w);
    end
    
    theta_vec = -varphi - 2*arg_gnT_w;
    
    P_recv = abs(f'*diag(exp(1j*theta_vec))*G*w)^2;
    sigma_noise = BS_conf.sigma_noise;      % Thermal noise.
    P_noise = sigma_noise^2;
    Rate = log2(1+P_recv/P_noise);  % In fact, this is the spectral efficiency.
end