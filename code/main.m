%% Construct MISO Single-user RIS Channel Estimation-Beamforming platform.
clear; clc;  
c = 299792458;
fc = 3.5e9;
lambda = c/fc;

Ny = 20;
Nz = 10;
N = Ny * Nz;
M = 8;              % Number of antennas at BS.

channel_type = 'rayleigh'; 
BL = 1000;              % Transmit blocklength.
Ep = 50;                % 50 unit energy for transmitting pilots. 
Ed = BL-Ep;             % energy for transmitting data. 
% (Ep) to joules: related to the transmitted power. 

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
Pt_UE = 100e-3;                 % 300mW UE transmit power.
sigma_noise = sqrt(P_noise);
RIS_conf.sigma_v = sqrt(PSD_noise * 100e6)*db2mag(10);     % 100M BW for IRF-sensing elements.

fprintf('Receiver noise power for each subcarrier = %.3f dBm\n', mag2db(sigma_noise)+30);
fprintf('RF sensor noise power = %.3f dBm\n', mag2db(RIS_conf.sigma_v)+30);
fprintf('User transmit power = %.3f dBm\n', pow2db(Pt_UE)+30);
fprintf('Typical path loss@200m = %.3f dB\n', mag2db(lambda/(4*pi*200)));
fprintf('Typical double-fading loss@(100m, 100m) = %.3f dB\n', mag2db((lambda/(4*pi*100))^2*N));

N_sim = 200;

Pt_BS_range     = flip(logspace(-1, 1, 10));
scan_len        = length(Pt_BS_range);
SE              = zeros(scan_len, 8);

fprintf('Start simulating...\n');
parfor idx_scan = 1:scan_len
    % Pt_BS = 1;      % Total transmit power is 1W. (distributed on M transmit antennas).
    Pt_BS = Pt_BS_range(idx_scan);
    
    L1 = 4;         % Number of BS-RIS NLOS paths.
    L2 = 4;         % Number of RIS-UE NLOS paths.
    kappa = 2;      % LOS ratio.
    Np_MMSE = N;         
    
    BS_conf         = struct();
    BS_conf.My      = M;            % Assume the BS is equipped with lambda/2 ULA.
    BS_conf.Mz      = 1;
    BS_conf.M       = BS_conf.My*BS_conf.Mz;
    BS_conf.Pt_BS   = Pt_BS;
    BS_conf.Pt_UE   = Pt_UE;
    BS_conf.BL      = BL;
    BS_conf.Ep      = Ep;
    BS_conf.Ed      = Ed;
    BS_conf.sigma_noise = sigma_noise;

    user_conf = struct();

    %% Run simulations.
    Rates = zeros(N_sim, 8);
    rng(0);
    for idx_sim =  1:N_sim
        [R_BS, theta_BS, psi_BS] = randPos(R_BS_range, theta_BS_range, psi_BS_range);
        pos_BS = R_BS*[cos(psi_BS)*cos(theta_BS), cos(psi_BS)*sin(theta_BS), sin(psi_BS)];
        BS_conf.pos_BS = pos_BS;    % For IRF methods, the algorithm needs to know where the BS is.
        BS_conf.BS_pos_sph = [R_BS, theta_BS, psi_BS];

        [R_UE, theta_UE, psi_UE] = randPos(R_UE_range, theta_UE_range, psi_UE_range);
        pos_UE = R_UE*[cos(psi_UE)*cos(theta_UE), cos(psi_UE)*sin(theta_UE), sin(psi_UE)];
        user_conf.pos_UE = pos_UE;
        user_conf.user_pos_sph = [R_UE, theta_UE, psi_UE]; 
        
        % Calculate the channel G, h(k), and the parameters alpha and beta.
        % (Based on the Friis transmission formula.)
        % The "channels" are complex power-based transfer functions.
        % Notice: The channel G must be highly structured and predictable.
        
        if strcmp(channel_type, 'SV')
            [f_LOS, G_LOS] = generate_channel_los(RIS_conf, BS_conf, [R_BS, theta_BS, psi_BS], [R_UE, theta_UE, psi_UE]);
            [f_NLOS, G_NLOS] = generate_channel_multipath([Ny, Nz], [M, 1], RIS_conf, pos_BS, [0,0,0], pos_UE, L1, L2);
            f = sqrt(kappa/(1+kappa))*f_LOS + sqrt(1/(1+kappa))*f_NLOS;
            G = sqrt(kappa/(1+kappa))*G_LOS + sqrt(1/(1+kappa))*G_NLOS;
        elseif strcmp(channel_type, 'rayleigh')
            [f, G] = generate_channel_Rayley(RIS_conf, BS_conf, user_conf);
        else
            error('channel_type not correct.');
        end

        % Baseline 1: Random phasing, random precoding...
        w = (randn([M,1])+1j*randn([M,1]))/sqrt(2);
        w = sqrt(Pt_BS)*w/norm(w);
        theta0 = exp(1j*2*pi*rand([N,1]));
        y = f'*diag(theta0)*G*w;
        Rate_random = log2(1+abs(y)^2/P_noise); % bit per transmitted symbol. 
        Rates(idx_sim, 1) = Rate_random;    

        % Baseline 2: Perform channel estimation by traditional methods: Orthogonal
        % pilots, LS-CE/LMMSE. (1) LMMSE for H;  
        [Y, Thetas] = transmit_pilot(RIS_conf, BS_conf, f, G, N); 

        [HT_hat, nmse_MF]   = estimate_HT(RIS_conf, BS_conf, Y, Thetas, f, G, N, 'MF');
        [w, theta]          = RIS_precode(RIS_conf, HT_hat, theta0);
        Rate_MF             = calc_rate(BS_conf, G, f, w, theta, N);
        Rates(idx_sim, 2)   = Rate_MF;  

        [HT_hat, nmse_LS]   = estimate_HT(RIS_conf, BS_conf, Y, Thetas, f, G, N, 'LS');
        [w, theta]          = RIS_precode(RIS_conf, HT_hat, theta0);
        Rate_LS             = calc_rate(BS_conf, G, f, w, theta, N);
        Rates(idx_sim, 3)   = Rate_LS;  
        
        [HT_hat, nmse_MMSE] = estimate_HT(RIS_conf, BS_conf, Y, Thetas, f, G, N, 'MMSE');
        [w, theta]          = RIS_precode(RIS_conf, HT_hat, theta0);
        Rate_MMSE           = calc_rate(BS_conf, G, f, w, theta, N);
        Rates(idx_sim, 4)   = Rate_MMSE;  

        % (2) MMSE for f only. Assume that G is known. 
        [f_hat, Np_LMMSE_f, ~]  = LMMSE_estimate_f(RIS_conf, BS_conf, f, G, Ep, Ed);
        [w, theta]              = RIS_precode(RIS_conf, G.'*diag(conj(f_hat)));     % H=diag(f*)G, here HT required. 
        Rate_LMMSE_f            = calc_rate(BS_conf, G, f, w, theta, Np_LMMSE_f);
        Rates(idx_sim, 5)       = Rate_LMMSE_f;
        
        % Use IRF to directly perform beamforming, using only 3 pilots.
        % Assume that G is known. 
        [P_recv_IRF, Rate_IRF] = IRF_CE_BF(RIS_conf, BS_conf, f, G, Ep, Ed, 3, channel_type);
        Rates(idx_sim, 6) = Rate_IRF;

        % Genie told me the G and f-channel + iterate between \theta and w. 
        [w, theta] = RIS_precode(RIS_conf, G.'*diag(conj(f)));
        Rates(idx_sim, 7) = calc_rate(BS_conf, G, f, w, theta, 0); 
    end
    
    SE(idx_scan, :) = mean(Rates);
    fprintf('Pt_{BS} = %.3f dB sim complete. %d/%d\n', pow2db(Pt_BS), idx_scan, scan_len);
end


%% Save Data. 
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

plot(dbP, SE(:,7), 'ko-.','MarkerSize',6, 'LineWidth', 2);  % Oracle
plot(dbP, SE(:,6), 'ro-','MarkerSize',6);                   % IRF
plot(dbP, SE(:,5), 'color', [0, 0, 1], 'LineStyle', '-', 'marker', 'x','MarkerSize',8);     % MMSE f
plot(dbP, SE(:,4), 'color', [1, 0, 1], 'LineStyle', '--', 'marker', 'p','MarkerSize',6);    % LMMSE H
plot(dbP, SE(:,2), 'color', [0, 1, 1], 'LineStyle', '-', 'marker', 's','MarkerSize',6);     % MF H
plot(dbP, SE(:,1), 'color', [1, 0, 0.9], 'LineStyle', '-', 'marker', 'x', 'MarkerSize',6);  % Random


set(gca,'FontName','Times New Roman');
grid on; box on;
legend('Oracle', 'Proposed-IRF + VM-EM', 'MMSE f', 'LMMSE H', 'MF H', 'Random');
xlabel('BS transmit power $P_{\rm max}$ (dBW)', 'interpreter', 'latex');
ylabel('$C_{\rm SI}$ (bps/Hz)', 'interpreter', 'latex');


%% Utilities.
function [R, theta, psi] = randPos(R_range, theta_range, psi_range)
    R       = R_range(1) + (R_range(2) - R_range(1))*rand();
    theta   = theta_range(1) + (theta_range(2) - theta_range(1))*rand();
    psi     = psi_range(1) + (psi_range(2) - psi_range(1))*rand();
end

