function [outputArg1,outputArg2] = CE_sparsity(N, M, K, f, G)
N1 = N(1); N2 = N(2); N = N1*N2;
M1 = M(1); M2 = M(2); M = M1*M2;
% unitary matrices
UN1=(1/sqrt(N1))*exp(-1i*2*pi*(0:N1-1)'*(0:N1-1)*(2/N1)/2);
UN2=(1/sqrt(N2))*exp(-1i*2*pi*(0:N2-1)'*(0:N2-1)*(2/N2)/2);
UN=kron(UN1,UN2);

UM1=(1/sqrt(M1))*exp(-1i*2*pi*(0:M1-1)'*(-(M1-1)/2:1:(M1/2))*(2/M1)/2);
UM2=(1/sqrt(M2))*exp(-1i*2*pi*(0:M2-1)'*(-(M2-1)/2:1:(M2/2))*(2/M2)/2);
UM=kron(UM1,UM2);

sample = 20;
Q_all = 32:16:128;
for s = 1:sample
    fprintf("s = %2d\n", s);
    H = zeros(N, M, K);
    for k = 1:K
        fk = f(:, k);
        Hk = G'*diag(fk);
        H(:,:,k)=(UM'*Hk*(UN.')')';
        %imagesc(abs(H(:,:,k)));
    end
    for i = 1:length(Q_all)
        Q = Q_all(i);
        Y = zeros(Q, M, K);
        W = ((rand(N,Q)>0.5)*2-1)/sqrt(N);
        A = (UN*W)';
        for k = 1:K
            noise = sqrt(sigma2)*(randn(Q,M)+1i*randn(Q,M))/sqrt(2);
            Y(:,:,k) = A*H(:,:,k)+noise;
        end
    end
end
end