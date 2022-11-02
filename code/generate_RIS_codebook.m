function Thetas = generate_RIS_codebook(N, Np)
    F = dftmtx(N);
    
    if Np <= N
        Thetas = F(:, 1:Np);    
    else
        Thetas = zeros(N, Np);
        Thetas(:,1:N) = F;
        Thetas(:, N+1:end) = exp(1j*2*pi*rand([N, Np-N]));
    end
end

