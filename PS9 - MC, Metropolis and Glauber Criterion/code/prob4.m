%% Problem 4 part (ii)
clear all;
close all;
clc;

d_list = 2:20;
f_list = zeros(1, length(d_list));
k_list = zeros(1, length(d_list));
for i = 1:length(d_list)
    echodemo('ho_glauber.m', 1);

    nsteps=100000;
    d = d_list(i);

    echodemo('ho_glauber.m', 2);

    nmax = 10;
    count = zeros(1, nmax);
    corr = zeros(1, nmax);
    for m = 1:nsteps-nmax
        for n = 1:nmax
            count(n) = count(n) + 1;
            corr(n) = corr(n) + xtraj(m) * xtraj(m + n);
        end
    end

    corr = corr ./ count;
    f_list(i) = num_accepted / h.count;
    lin_fit_params = polyfit(1:nmax, log(abs(corr)), 1);
    k_list(i) = -lin_fit_params(1);
    
    subplot(4, 5, i);
    plot(1:nmax, log(abs(corr)));
    title(sprintf('d = %.2f', d));
    ylim([-4, 0]);
end


figure;
subplot(2, 1, 1)
plot(d_list, 1./k_list);
title('n_{corr} vs. d')
xlabel('d')
ylabel('n_{corr}')

subplot(2, 1, 2)
plot(f_list, 1./k_list)
title('n_{corr} vs. f_{acc}')
xlabel('f_{acc}')
ylabel('n_{corr}')