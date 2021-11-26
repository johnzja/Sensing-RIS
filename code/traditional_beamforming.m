function [w, theta] = traditional_beamforming(f, G, P)
% P - transmit power
% f - N*1
% G - N*M
% w - M*1
% theta - N*1
[N, ~] = size(G);
% objective: max |f^H*Theta*G*w|^2
threshold = 0.01;
iter = 10;
objective = zeros(2*iter, 1);
theta = exp(1j*2*pi*rand(N, 1));
for i = 1:iter
    % update w
    w = G'*diag(theta)'*f;
    w = sqrt(P)*w/norm(w);
    objective(2*i-1) = abs(f'*diag(theta)*G*w)^2;
    % update theta
    theta = exp(-1j*angle(diag(conj(f))*G*w));
    objective(2*i) = abs(f'*diag(theta)*G*w)^2;
    if abs(objective(2*i)-objective(2*i-1))/objective(2*i) < threshold
        break;
    end
end
end

