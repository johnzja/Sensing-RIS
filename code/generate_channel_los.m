function [f, G] = generate_channel_los(N, M, BS_pos, RIS_pos, user_pos, fc)
N = N(1)*N(2);
M = M(1)*M(2);
% distance calculation
f_dis = norm(RIS_pos-user_pos);
G_dis = norm(BS_pos-RIS_pos);
% f - N*1
% G - N*M
lambda = 299792458/fc;
% large-scale fading
f_fading = lambda/(4*pi)/f_dis; % Friis Transmission Formula.
G_fading = lambda/(4*pi)/G_dis;
% small-scale fading
f = (randn(N, 1)+1j*randn(N, 1))/sqrt(2);
G = (randn(N, M)+1j*randn(N, M))/sqrt(2);

f = f*f_fading;
G = G*G_fading;
end

