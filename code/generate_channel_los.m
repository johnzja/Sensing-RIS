function [f, G] = generate_channel_los(RIS_conf, BS_conf, BS_pos_sph, user_pos_sph)
    Ny = RIS_conf.Ny; Nz = RIS_conf.Nz;
    My = BS_conf.My; Mz = BS_conf.Mz;

    % distance calculation. Assume that RIS is located at [0,0,0].
    f_dis = user_pos_sph(1);   
    G_dis = BS_pos_sph(1);

    % f - N*1
    % G - N*M
    lambda = RIS_conf.lambda;
    % large-scale fading
    f_fading = lambda/(4*pi)/f_dis;        % Friis Transmission Formula. (1, K).
    G_fading = lambda/(4*pi)/G_dis;

    % Formula: exp(-j*(nz*sin(psi)+ny*cos(psi)*sin(theta)) for RIS.
    % Assume that the BS is always faced to RIS.

    psi_BS = BS_pos_sph(3);
    theta_BS = BS_pos_sph(2);

    phase_G = kron(exp(1j*pi*sin(psi_BS)*(0:Nz-1)), exp(1j*pi*cos(psi_BS)*sin(theta_BS)*(0:Ny-1)));
    phase_BS = kron(exp(-1j*pi*sin(psi_BS)*(0:Mz-1)), exp(-1j*pi*cos(psi_BS)*sin(theta_BS)*(0:My-1)));
    % f  = (randn(N, 1)+1j*randn(N, 1))/sqrt(2)*f_fading;
    % G = (randn(N, M)+1j*randn(N, M))/sqrt(2);
    % G = G*G_fading;
    G = kron(phase_G.', phase_BS);  % size(G) = [N, M].
    G = G*G_fading;                 % Large-scale fading.

    theta_UE = user_pos_sph(2);
    psi_UE = user_pos_sph(3);
    phase_f = kron(exp(1j*pi*sin(psi_UE)*(0:Nz-1)), exp(1j*pi*cos(psi_UE)*sin(theta_UE)*(0:Ny-1)));
    f = phase_f.' * f_fading;
end

