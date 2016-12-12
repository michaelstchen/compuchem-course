%% Problem 1 parts (ii) + (iii)
clear all;
close all;
clc;

Nside=6;
density = 0.7;

echodemo('hard_disks.m', 2);
echodemo('hard_disks.m', 3)

figure
plot(h.vals,h.hist./(nsweeps*pi*N*density*h.binwidth*h.vals))
xlabel('r')
ylabel('g(r)')

%% Problem 1 part (iv)
clear all;
close all;
clc;

dlist = 0.1:0.1:0.8;
glist = zeros(1, length(dlist));
for di = 1:length(dlist)
    clearvars -except dlist glist di
    density = dlist(di);
    Nside = 6;
    
    echodemo('hard_disks.m', 2);
    echodemo('hard_disks.m', 3);
    
    subplot(3, 3, di);
    gr = h.hist./(nsweeps*pi*N*density*h.binwidth*h.vals);
    plot(h.vals,gr)
    title(sprintf('density = %.1f', density))
    xlabel('r')
    ylabel('g(r)')
    xlim([0 5])
    
    glist(di) = gr(11);
end

fprintf('Density & g(sigma+) \\\\\n')
for ig = 1:length(glist)
    fprintf('%.1f & %.3f \\\\\n', dlist(ig), glist(ig))
end

%% Problem 1 part (v)

p = @(g, rho) (rho + pi * rho^2 * g / 2);
plist = zeros(1, length(glist));
for pi = 1:length(plist)
    plist(pi) = p(glist(pi), dlist(pi));
end

plot(dlist, plist);
xlabel('density')
ylabel('pressure')