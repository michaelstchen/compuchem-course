%% Problem 1 part (iii)
clear all;
close all;
clc;

h = struct('numbins', 0,...
           'range', [1,10],...
           'hist', [],...
           'vals', [],...
           'count', 0);
       
err = @(prob, n_total) sqrt((prob .* (1 - prob)) / n_total);
       
for i = 1:15000
    h = histo_int(h, randint(1, 10));
end

prob_dist = h.hist / h.count;

hold on;
bar(h.vals, prob_dist);
title('RandInt Distrib For 15000 samples');
ylabel('Sampled Probability');
xlabel('Integer');
% errorbar(h.vals, prob_dist, err(prob_dist, h.count), 'o');
hold off;


%% Problem 1 part (iv)
clear all;
close all;
clc;

h = struct('numbins', 0,...
           'range', [1,10],...
           'hist', [],...
           'vals', [],...
           'count', 0);
       
err = @(prob, n_total) sqrt((prob .* (1 - prob)) / n_total);

n_samps = [1000000];
for i = 1:length(n_samps)
    n_total = n_samps(i);
    h = struct('numbins', 0,...
           'range', [1,10],...
           'hist', [],...
           'vals', [],...
           'count', 0);
       
    for j = 1:n_total
        h = histo_int(h, randint(1, 10));
    end
    
    prob_dists = h.hist / h.count;
    
    figure;
    hold on;
    bar(h.vals, prob_dists);
    title(sprintf('RandInt Distrib For %d samples', n_total));
    ylabel('Sampled Probability');
    xlabel('Integer');
    errorbar(h.vals, prob_dists, err(ones(1, 10)/10, h.count), 'o');
    hold off;
end