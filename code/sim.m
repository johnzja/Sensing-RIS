%% Construct MISO Single-user RIS Channel Estimation-Beamforming platform.
clear; clc;
c = 299792458;
fc = 10e9;
lambda = c/fc;

Ny = 10;
Nz = 10;
N = Ny * Nz;
M = 10;             % Number of antennas at BS.

% lambda/2-spaced RIS.
RIS_conf.Ny = Ny;
RIS_conf.Nz = Nz;
RIS_conf.N = N;
RIS_conf.center_y = ((1:Ny)-mean((1:Ny)))*lambda/2;
RIS_conf.center_z = ((1:Nz)-mean((1:Nz)))*lambda/2;
RIS_conf.lambda = lambda;

% Setup hyperparameters for the location of BS and Users.
R_BS_range = [20, 100];
theta_BS_range = [-pi/2, pi/2];
psi_BS_range = [-pi/3, pi/3];

R_UE_range = [10, 200];
theta_UE_range = [-pi/2, pi/2];
psi_UE_range = [-pi/3, pi/3];

% Setup simulation necessities.
N_sim = 1000;
rng(0);
Pt_BS = 1;      % Total transmit power is 1W. (distributed on M transmit antennas).
Pt_UE = 300e-3; % 300mW UE transmit power.
L1 = 4;         % Number of BS-RIS NLOS paths.
L2 = 4;         % Number of RIS-UE NLOS paths.
kappa = 0.8;    % LOS ratio.
Np = 10;        % Number of pilots in traditional beamforming.

BS_conf.M = M;
BS_conf.Pt_BS = Pt_BS;


%% Run simulation.
for idx_sim =  1:N_sim
    [R_BS, theta_BS, psi_BS] = randPos(R_BS_range, theta_BS_range, psi_BS_range);
    pos_BS = R_BS*[cos(psi_BS)*cos(theta_BS), cos(psi_BS)*sin(theta_BS), sin(psi_BS)];
    BS_conf.pos_BS = pos_BS;    % For IRF methods, the algorithm needs to know where the BS is.
    BS_conf.BS_pos_sph = [R_BS, theta_BS, psi_BS];
    pos_UE = zeros(K, 3);

    [R_UE, theta_UE, psi_UE] = randPos(R_UE_range, theta_UE_range, psi_UE_range);
    pos_UE(k, :) = R_UE*[cos(psi_UE)*cos(theta_UE), cos(psi_UE)*sin(theta_UE), sin(psi_UE)];

    % Calculate the channel G, h(k), and the parameters alpha and beta.
    % (Based on the Friis transmission formula.)
    % The "channels" are complex power-based transfer functions.
    % Notice: The channel G must be highly structured and predictable.
    [f_LOS, G_LOS] = generate_channel_los([Ny, Nz], [M, 1], RIS_conf, [R_BS, theta_BS, psi_BS], [R_UE, theta_UE, psi_UE]);
    [f_NLOS, G_NLOS] = generate_channel_multipath([Ny, Nz], [M, 1], RIS_conf, pos_BS, [0,0,0], pos_UE, L1, L2);
    f = sqrt(kappa/(1+kappa))*f_LOS + sqrt(1/(1+kappa))*f_NLOS;
    G = sqrt(kappa/(1+kappa))*G_LOS + sqrt(1/(1+kappa))*G_NLOS;
    
    % Perform channel estimation by traditional methods: Orthogonal
    % pilots, MMSE. Assume the RIS to be continuously adjustable within
    % [0,2pi].
    [P_recv_traditional, Rate_traditional] = ...
        traditional_CE_BF(N, M, RIS_conf, BS_conf, f, G, Np);

    % TODO: Use IRF to directly perform beamforming, using the same amount
    % of pilots.
    
end


%% Utilities.
function [R, theta, psi] = randPos(R_range, theta_range, psi_range)
    R       = R_range(1) + (R_range(2) - R_range(1))*rand();
    theta   = theta_range(1) + (theta_range(2) - theta_range(1))*rand();
    psi     = psi_range(1) + (psi_range(2) - psi_range(1))*rand();
end

