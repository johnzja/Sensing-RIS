%% Setup parameters for CRLB. 
L = 256;

gamma_bar_list = 0:5:100;
K_list = 0:0.05:1;
psi_arr = cos((0:(L-1))/L*2*pi);
sin_psi_arr_squared = sin(psi_arr).^2;
cos_psi_arr = cos(psi_arr);

recip_CRLB = zeros(length(gamma_bar_list), length(K_list));

g = @(x)(0.25*sqrt(pi./x).*exp(-x/2).*((1+1./x).*besseli(0,x/2)+besseli(1,x/2)));
relu = @(x)(max(zeros(size(x)), x));

for idx_gamma = 1:length(gamma_bar_list)
    for idx_K = 1:length(K_list)
        gamma_bar = gamma_bar_list(idx_gamma);
        K = K_list(idx_K);
        gamma = gamma_bar*(1+K*cos_psi_arr);    % calculate the real gamma_l. 
        % tmp = K^2 * gamma_bar^2 * sum(sin_psi_arr_squared .* relu(1./gamma - g(gamma)));
        tmp = K^2 * gamma_bar^2 * sum(sin_psi_arr_squared .* get_h(gamma));
        recip_CRLB(idx_gamma, idx_K) = tmp;
    end
end

%% Plot reciprocal CRLB.
set(0,'DefaultLineMarkerSize',4);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth',1.4);
set(0,'defaultfigurecolor','w');

[X,Y] = meshgrid(gamma_bar_list, K_list);
surf(X', Y', recip_CRLB);
view(0,90);
colorbar;

set(gca,'FontName','Times New Roman');
xlabel('$\bar{\gamma}$', 'interpreter', 'latex')
ylabel('$K$', 'interpreter', 'latex')
zlabel('Precision $1/{\rm CRLB}(\varphi)$ (dB)', 'interpreter', 'latex')


%% Utils.
function h=get_h(gamma)
    persistent g_values;
    persistent gamma_range;
    % Load the data.
    if ~exist('g_values', 'var') || isempty(g_values)
        load data/preciseCRLB_data.mat Expression1;
        gamma_range = 0.1:0.01:5;
        g_values = Expression1;
    end

    tmp = zeros(size(gamma));
    for idx = 1:length(gamma)
        if 0<=gamma(idx) && gamma(idx) <= gamma_range(end)
            tmp(idx) = interp1(gamma_range, g_values, gamma(idx), 'spline', 'extrap');
        else
            tmp(idx) = 1/gamma(idx) - sqrt(pi/gamma(idx)) * exp(-gamma(idx)/2) * ((1+1/gamma(idx))*besseli(0,gamma(idx)/2) + besseli(1,gamma(idx)/2))/4;
        end
    end
    h=tmp;
end