function [f_hat, G_hat] = traditional_channel_estimation(N, M, RIS_conf, f, G, Np)
    Ny = N(1); Nz = N(2);
    My = M(1); Mz = M(2);
    N = Ny * Nz;
    
    F = dftmtx(N);
    Thetas = F(:, 1:Np);
    
    % Channel estimation scheme: Just estimate fk and G for each user.


end