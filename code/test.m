clear all;
N = 32; % RIS element
M = 4; % BS antenna
BS_pos = [0, 0];
RIS_pos = [5, 5];
user_pos = [50, 0];
fc = 5e8;
[f, G] = generate_channel(N, M, BS_pos, RIS_pos, user_pos, fc);
P = 1;
[w, theta] = traditional_beamforming(f, G, P);