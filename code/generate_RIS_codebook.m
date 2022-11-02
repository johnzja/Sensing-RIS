function Thetas = generate_RIS_codebook(N, Np, alg_mode)
    switch alg_mode
        case 'DFT'
            F = dftmtx(N);
            
            if Np <= N
                Thetas = F(:, 1:Np);    
            else
                Thetas = zeros(N, Np);
                Thetas(:,1:N) = F;
                Thetas(:, N+1:end) = exp(1j*2*pi*rand([N, Np-N]));
            end
        case 'random'
            Thetas = exp(1j*2*pi*rand([N, Np]));
        otherwise
            error('unknown alg_mode'); 
    end
end

