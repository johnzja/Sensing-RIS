%% Construct MISO-Multiuser RIS Channel Estimation-Beamforming platform.

c = 299792458;
fc = 10e9;
lambda = c/fc;

Ny = 10;
Nz = 10;
N = Ny * Nz;
K = 10;             % Number of UEs.
M = 10;             % Number of antennas at BS.

% lambda/2-spaced RIS.
RIS_conf.Ny = Ny;
RIS_conf.Nz = Nz;
RIS_conf.N = N;
RIS_conf.center_y = ((1:Ny)-mean((1:Ny)))*lambda/2;
RIS_conf.center_z = ((1:Nz)-mean((1:Nz)))*lambda/2;

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
Pt_BS = 1;      % Total transmit power is 1W. (on M transmit antennas).
Pt_UE = 300e-3; % 300mW UE transmit power.

R_UE = zeros(K, 1);
theta_UE = zeros(K, 1);
psi_UE = zeros(K, 1);

%% Run simulation.
for idx_sim =  1:N_sim
    [R_BS, theta_BS, psi_BS] = randPos(R_BS_range, theta_BS_range, psi_BS_range);
    for k=1:K
        [r_ue, theta_ue, psi_ue] = randPos(R_UE_range, theta_UE_range, psi_UE_range);
        R_UE(k) = r_ue; theta_UE(k) = theta_ue; psi_UE(k) = psi_ue;
    end
    
    % Calculate the channel G, h(k), and the parameters alpha and beta.
    % Based on the Friis transmission formula.
    % The "channels" are complex power-based transfer functions.
    G = zeros(M, N);
    for m = 1:M
        for n = 1:N
            nx = mod(n-1, Nx) + 1;
            nz = floor((n-1)/Nx)+1;
            G(m,n) = (lambda/(4*pi*
        end
    end
    
end


%% Utilities.
function [R, theta, psi] = randPos(R_range, theta_range, psi_range)
    R       = R_range(1) + (R_range(2) - R_range(1))*rand();
    theta   = theta_range(1) + (theta_range(2) - theta_range(1))*rand();
    psi     = psi_range(1) + (psi_range(2) - psi_range(1))*rand();
end

