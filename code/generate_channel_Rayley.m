function [f, G] = generate_channel_Rayley(RIS_conf, BS_conf, user_conf)
    N = RIS_conf.N;
    M = BS_conf.M;

    % distance calculation. Assume that RIS is located at [0,0,0].
    f_dis = user_conf.user_pos_sph(1);   
    G_dis = BS_conf.BS_pos_sph(1);

    % f - N*1
    % G - N*M
    lambda = RIS_conf.lambda;
    % large-scale fading
    f_fading = lambda/(4*pi)/(f_dis^(2.0/2));
    G_fading = lambda/(4*pi)/(G_dis^(2.0/2));

    G = (randn([N, M]) + 1j*randn([N, M])) /sqrt(2)*G_fading;
    f = (randn(N, 1) + 1j*randn(N, 1)) /sqrt(2)*f_fading;
end