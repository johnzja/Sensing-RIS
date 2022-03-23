%% Construct MISO Single-user RIS Channel Estimation-Beamforming platform.
clear; clc;
c = 299792458;
fc = 10e9;
lambda = c/fc;

Ny = 40;
Nz = 30;
N = Ny * Nz;
M = 8;              % Number of antennas at BS.

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
PSD_noise = db2pow(-174-30);    % Setup the PSD of the thermo-noise.
BW = 180e3;                     % System baseband BW on each subcarrier = 180kHz.
P_noise = PSD_noise * BW;       % Thermo-noise for BS receiver equipped with M RF-chains.
Pt_UE = 300e-3;                 % 300mW UE transmit power.
sigma_noise = sqrt(P_noise);
RIS_conf.sigma_v = sqrt(PSD_noise * 100e6);     % 100M BW for IRF-sensing elements.
    
N_sim = 200;

Pt_BS_range = logspace(-1,1, 10);
SE = zeros(length(Pt_BS_range), 4);

for idx_scan = 1:length(Pt_BS_range)
    % Pt_BS = 1;      % Total transmit power is 1W. (distributed on M transmit antennas).
    Pt_BS = Pt_BS_range(idx_scan);
    
    L1 = 4;         % Number of BS-RIS NLOS paths.
    L2 = 4;         % Number of RIS-UE NLOS paths.
    kappa = 2;    % LOS ratio.
    Np = N;         % Number of pilots in traditional beamforming.
    
    BS_conf         = struct();
    BS_conf.My      = M;            % Assume the BS is equipped with lambda/2 ULA.
    BS_conf.Mz      = 1;
    BS_conf.M       = BS_conf.My*BS_conf.Mz;
    BS_conf.Pt_BS   = Pt_BS;
    BS_conf.Pt_UE   = Pt_UE;
    BS_conf.sigma_noise = sigma_noise;

    %% Run simulations.
    Rates = zeros(N_sim, 4);
    rng(0);
    for idx_sim =  1:N_sim
        [R_BS, theta_BS, psi_BS] = randPos(R_BS_range, theta_BS_range, psi_BS_range);
        pos_BS = R_BS*[cos(psi_BS)*cos(theta_BS), cos(psi_BS)*sin(theta_BS), sin(psi_BS)];
        BS_conf.pos_BS = pos_BS;    % For IRF methods, the algorithm needs to know where the BS is.
        BS_conf.BS_pos_sph = [R_BS, theta_BS, psi_BS];

        [R_UE, theta_UE, psi_UE] = randPos(R_UE_range, theta_UE_range, psi_UE_range);
        pos_UE = R_UE*[cos(psi_UE)*cos(theta_UE), cos(psi_UE)*sin(theta_UE), sin(psi_UE)];

        % Calculate the channel G, h(k), and the parameters alpha and beta.
        % (Based on the Friis transmission formula.)
        % The "channels" are complex power-based transfer functions.
        % Notice: The channel G must be highly structured and predictable.
        [f_LOS, G_LOS] = generate_channel_los(RIS_conf, BS_conf, [R_BS, theta_BS, psi_BS], [R_UE, theta_UE, psi_UE]);
        [f_NLOS, G_NLOS] = generate_channel_multipath([Ny, Nz], [M, 1], RIS_conf, pos_BS, [0,0,0], pos_UE, L1, L2);
        f = sqrt(kappa/(1+kappa))*f_LOS + sqrt(1/(1+kappa))*f_NLOS;
        G = sqrt(kappa/(1+kappa))*G_LOS + sqrt(1/(1+kappa))*G_NLOS;

        % Baseline 1: Random phasing, random precoding...
        w = (randn([M,1])+1j*randn([M,1]))/sqrt(2);
        w = sqrt(Pt_BS)*w/norm(w);
        theta = exp(1j*2*pi*rand([N,1]));
        y = f'*diag(theta)*G*w;
        Rate_random = log2(1+abs(y)^2/P_noise);
        Rates(idx_sim, 1) = Rate_random;

        % Baseline 2: Perform channel estimation by traditional methods: Orthogonal
        % pilots, LS-CE/MMSE. Assume the RIS to be continuously adjustable within
        % [0,2pi].
        [P_recv_traditional, Rate_traditional] = traditional_CE_BF(RIS_conf, BS_conf, f, G, Np);
        Rates(idx_sim, 2) = Rate_traditional;

        % Use IRF to directly perform beamforming, using only 3 pilots.
        [P_recv_IRF, Rate_IRF] = IRF_CE_BF(RIS_conf, BS_conf, f, G, 3);
        Rates(idx_sim, 3) = Rate_IRF;

        % Genie told me the channel.
        [P_recv_Oracle, Rate_Oracle] = traditional_BF(RIS_conf, BS_conf, f, G);
        Rates(idx_sim, 4) = Rate_Oracle;
    end
    R = mean(Rates);
    SE(idx_scan, :) = R;
    fprintf('Pt_{BS} = %f \n', Pt_BS);
end
disp('sim complete. Files saved at data/.');
save('data/IRF_sim.mat', 'Pt_BS_range', 'SE', 'M', 'Pt_UE', 'RIS_conf', 'sigma_noise', 'N_sim');


%% Analyze the data.
load('data/IRF_sim.mat');
dbP = 10*log10(Pt_BS_range).';    % in dBW.

set(0,'DefaultLineMarkerSize',6);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth',1.4);
set(0,'defaultfigurecolor','w');

figure('color',[1 1 1]); hold on;
plot(dbP, SE(:,1), 'color', [1, 0, 0.9], 'LineStyle', '-', 'marker', 'x');
plot(dbP, SE(:,2), 'color', [0, 0, 1], 'LineStyle', '-', 'marker', 'x');
plot(dbP, SE(:,3), 'ro-');
plot(dbP, SE(:,4), 'ko-.');

set(gca,'FontName','Times New Roman');
grid on; box on;
legend('Random', 'MMSE', 'Proposed-IRF', 'Oracle');
xlabel('BS transmit power (dBW)', 'interpreter', 'latex');
ylabel('Capacity (bps/Hz)', 'interpreter', 'latex');


%% Utilities.
function [R, theta, psi] = randPos(R_range, theta_range, psi_range)
    R       = R_range(1) + (R_range(2) - R_range(1))*rand();
    theta   = theta_range(1) + (theta_range(2) - theta_range(1))*rand();
    psi     = psi_range(1) + (psi_range(2) - psi_range(1))*rand();
end

