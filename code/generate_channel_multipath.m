function [f, G] = generate_channel_multipath(N, M, RIS_conf, BS_pos, RIS_pos, user_pos, L1, L2)
N1 = N(1); N2 = N(2); N = N1*N2; sz_N = [N2, N1];
M1 = M(1); M2 = M(2); M = M1*M2; sz_M = [M2, M1];
% distance calculation
f_dis = norm(RIS_pos-user_pos); 
G_dis = norm(BS_pos-RIS_pos);
% large-scale fading
C0 = RIS_conf.lambda/(4*pi);    % Friis-transmission formula, coefficient.
p_G = C0*G_dis^(-2.2);
p_f = C0*f_dis.^(-2.8);

% SV model
% L1: Number of paths from BS to RIS.
% L2: Number of paths from RIS to UE.

% Setup (logical) angles. Note: This angle distribution is not uniform in
% the physical domain.
M1index = ((1:M1) - (1+M1)/2)*(2/(M1));
M2index = ((1:M2) - (1+M2)/2)*(2/(M2));
N1index = ((1:N1) - (1+N1)/2)*(2/(N1));
N2index = ((1:N2) - (1+N2)/2)*(2/(N2));

% generate the physical angles of G
index = randperm(N); 
[y, x] = ind2sub(sz_N, index(1:L1));
phi1 = N1index(x);
phi2 = N2index(y);

index = randperm(M); 
[y, x] = ind2sub(sz_M, index(1:L1));
psi1 = M1index(x);
psi2 = M2index(y);

% generate the gains of G
alpha_G = (randn([L1, 1]) + 1j*randn([L1, 1])) / sqrt(2);

% generate G channel (N*M)
G = zeros(M, N);
for s = 1:L1
    a1 = 1/sqrt(N1)*exp(-1i*2*pi*(0:N1-1)'*phi1(s)/2);  % Here: /2 means lambda/2.
    a2 = 1/sqrt(N2)*exp(-1i*2*pi*(0:N2-1)'*phi2(s)/2);
    a = kron(a1,a2);
    b1 = 1/sqrt(M1)*exp(-1i*2*pi*(0:M1-1)'*psi1(s)/2);
    b2 = 1/sqrt(M2)*exp(-1i*2*pi*(0:M2-1)'*psi2(s)/2);
    b = kron(b1,b2);
    G = G + alpha_G(s)*b*a.';
end
G = sqrt(p_G*M*N/L1)*G';

% generate the physical angles of f.
% Note here: Xiuhong introduces "common paths" here. But it is unnecessary
% for our design.
% generate the channel for each user k

alpha = (randn([L2, 1]) + 1j*randn([L2, 1])) / sqrt(2);
fk=zeros(N,1);

index = randperm(N);
[y, x] = ind2sub(sz_N, index(1:L2));      
phi1(1:L2) = N1index(x);
phi2(1:L2) = N2index(y);
for s = 1:L2
    a1 = 1/sqrt(N1)*exp(-1i*2*pi*(0:N1-1)'*phi1(s)/2);
    a2 = 1/sqrt(N2)*exp(-1i*2*pi*(0:N2-1)'*phi2(s)/2);
    a = kron(a1,a2);
    fk = fk + alpha(s)*a;
end
f = sqrt(p_f(k)*N/L2)*fk; 

end

