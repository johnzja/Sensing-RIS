gamma_range = 0.1:0.01:5;
load data/preciseCRLB_data.mat;

g_values = Expression1;

plot(gamma_range, g_values); grid on;

%% 
x = [1,2,3,4, 5, 6];  % Assume the input x to be a vector.
h = zeros(size(x));
for idx = 1:length(x)
    if 0<=x(idx) && x(idx) <= gamma_range(end)
        h(idx) = interp1(gamma_range, g_values, x(idx), 'spline', 'extrap');
    else
        h(idx) = 1/x(idx) - sqrt(pi/x(idx)) * exp(-x(idx)/2) * ((1+1/x(idx))*besseli(0,x(idx)/2) + besseli(1,x(idx)/2))/4;
    end
end
disp(h);
