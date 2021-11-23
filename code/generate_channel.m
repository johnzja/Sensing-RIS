function [f, G] = generate_channel(N, M, BS_pos, RIS_pos, user_pos, fc)
% distance calculation
f_dis = norm(RIS_pos-user_pos);
G_dis = norm(BS_pos-RIS_pos);
% f - N*1
% G - N*M
lambda = 3e8/fc;
% large-scale fading
f_fading = sqrt(lambda^2/16/pi)/f_dis;
G_fading = sqrt(lambda^2/16/pi)/G_dis;
% small-scale fading
f = (randn(N, 1)+1j*randn(N, 1))/sqrt(2);
G = (randn(N, M)+1j*randn(N, M))/sqrt(2);

f = f*f_fading;
G = G*G_fading;
end

