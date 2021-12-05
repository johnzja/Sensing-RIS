clear all;
N = [8, 8]; % RIS element
M = [2, 2]; % BS antenna
K = 1;
BS_pos = [0, 0];
RIS_pos = [5, 5];
user_pos = [50, 0];
fc = 5e8;
L1 = 2; L2 = 8; Lc = 0;
[f, G] = generate_channel_multipath(N, M, K, BS_pos, RIS_pos, user_pos, L1, L2, Lc);
P = 1;
[w, theta] = traditional_beamforming(f, G, P);