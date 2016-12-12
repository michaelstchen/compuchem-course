%% Problem 1 part (ii)
clear all;
close all;
clc;

echodemo('ho_metro.m', 1);

nsteps=100000;
d = 1;

echodemo('ho_metro.m', 2);

hold on;
bar(h.vals,h.hist/(h.count*h.binwidth))
plot(h.vals,exp(-0.5*(h.vals .^2)/T)/sqrt(2*pi*T), '-r')
title('Histogram of Sampled Positions');
xlabel('x');
ylabel('normalized count');
hold off;

f_acc = num_accepted / h.count;


%% Problem 2 part (iii)
clear all;
close all;
clc;

d_list = [0.01, 0.5];
d_list = cat(2, d_list, 1:20);
f_list = zeros(1, length(d_list));
for i = 1:length(d_list)

    echodemo('ho_metro.m', 1);

    nsteps=50000;
    d = d_list(i);

    echodemo('ho_metro.m', 2);

    subplot(5, 5, i);
    hold on;
    plot(h.vals,h.hist/(h.count*h.binwidth))
    %plot(h.vals,exp(-0.5*(h.vals .^2)/T)/sqrt(2*pi*T))
    title(sprintf('d = %.2f', d));
    hold off;

    f_list(i) = num_accepted / h.count;

end

fprintf('d & f_acc \\\\ \n');
for i = 1:length(f_list)
    fprintf('%.2f & %.4f \\\\ \n', d_list(i), f_list(i));
end


%% Problem 2 part (v)
clear all;
close all;
clc;

% efficient sampling
echodemo('ho_metro.m', 1);

nsteps=10000;
d = 5;

echodemo('ho_metro.m', 2);

subplot(3, 1, 1)
plot(1:nsteps,xtraj);
title('Trajectory (Efficient d=5)');
xlabel('iteration');
ylabel('position');

% small d
echodemo('ho_metro.m', 1);

nsteps=10000;
d = 0.01;

echodemo('ho_metro.m', 2);

subplot(3, 1, 2)
plot(1:nsteps,xtraj);
title('Trajectory (Small d=0.01)');
xlabel('iteration');
ylabel('position');

% large d
echodemo('ho_metro.m', 1);

nsteps=10000;
d = 20;

echodemo('ho_metro.m', 2);

subplot(3, 1, 3)
plot(1:nsteps,xtraj);
title('Trajectory (Large d=20)');
xlabel('iteration');
ylabel('position');