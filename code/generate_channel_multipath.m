function [f, G] = generate_channel_multipath(N, M, K, BS_pos, RIS_pos, user_pos, L1, L2, Lc)
N1 = N(1); N2 = N(2); N = N1*N2;
M1 = M(1); M2 = M(2); M = M1*M2;
% distance calculation
f_dis = vecnorm((RIS_pos-user_pos)')'; % K*1 vector
G_dis = norm(BS_pos-RIS_pos);
% large-scale fading
C0 = 0.001;
p_G = C0*G_dis^(-2.2);
p_f = C0*f_dis.^(-2.8);
% SV model
M1index = (-(N1-1)/2:1:(M1/2))'*(2/M1);
M2index = (-(N2-1)/2:1:(M2/2))'*(2/M2);
N1index = (-(N1-1)/2:1:(N1/2))'*(2/N1);
N2index = (-(N2-1)/2:1:(N2/2))'*(2/N2);

% generate the physical angles of G
index = randperm(N); 
x = ceil(index(1:L1)/N2);
y = index(1:L1)-N2*(x-1);
phi1 = N1index(x);
phi2 = N2index(y);

index = randperm(M); 
x = ceil(index(1:L1)/M2);
y = index(1:L1)-M2*(x-1);
psi1 = M1index(x);
psi2 = M2index(y);

% generate the gains of G
alpha_G = (normrnd(0, 1, L1, 1) + 1i*normrnd(0, 1, L1, 1)) / sqrt(2);

% generate G channel (N*M)
G = zeros(M, N);
for l = 1:L1
    a1 = 1/sqrt(N1)*exp(-1i*2*pi*(0:N1-1)'*phi1(l)/2);
    a2 = 1/sqrt(N2)*exp(-1i*2*pi*(0:N2-1)'*phi2(l)/2);
    a = kron(a1,a2);
    b1 = 1/sqrt(M1)*exp(-1i*2*pi*(0:M1-1)'*psi1(l)/2);
    b2 = 1/sqrt(M2)*exp(-1i*2*pi*(0:M2-1)'*psi2(l)/2);
    b = kron(b1,b2);
    G = G + alpha_G(l)*b*a.';
end
G = sqrt(p_G*M*N/L1)*G';

% generate the physical angles of f
f = zeros(N,K);
index = randperm(N); 
x = ceil(index(1:Lc)/N2);
y = index(1:Lc)-N2*(x-1);
phi1c = N1index(x);
phi2c = N2index(y);

% generate the channel for each user k

for k = 1:K
    alpha = (normrnd(0, 1, L2, 1) + 1i*normrnd(0, 1, L2, 1)) / sqrt(2);
    fk=zeros(N,1);

    phi1(1:Lc) = phi1c;
    phi2(1:Lc) = phi2c;
    index = randperm(N);
    x = ceil(index(1:L2-Lc)/N2);
    y = index(1:L2-Lc)-N2*(x-1);       
    phi1(Lc+1:L2) = N1index(x);
    phi2(Lc+1:L2) = N2index(y);
    for l = 1:L2
        a1 = 1/sqrt(N1)*exp(-1i*2*pi*(0:N1-1)'*phi1(l)/2);
        a2 = 1/sqrt(N2)*exp(-1i*2*pi*(0:N2-1)'*phi2(l)/2);
        a = kron(a1,a2);
        fk = fk + alpha(l)*a;
    end
    f(:, k)= sqrt(p_f(k)*N/L2)*fk; 
end
end

