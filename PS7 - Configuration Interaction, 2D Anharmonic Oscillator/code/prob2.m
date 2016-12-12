%% part (v)
clear all;
close all;
clc;

fprintf('\nGround State Energy: %.4f\n', CIOscfunc(4, 1));

%% part (vi)
clear all;
close all;
clc;

a_list = 0:0.2:10;
N_orb_list = [2, 5, 10];
N_orb_energs = zeros(length(N_orb_list), length(a_list));
for N_orb_i = 1:length(N_orb_list)
    N_orb = N_orb_list(N_orb_i);
    for ai = 1:length(a_list)
        a = a_list(ai);
        N_orb_energs(N_orb_i, ai) = CIOscfunc(N_orb, a);
    end
end

% run('../PS5/prob3.m');
close all;

hold on;
plot(a_list, N_orb_energs(1,:), 'g');
plot(a_list, N_orb_energs(2,:), 'm');
plot(a_list, N_orb_energs(3,:), 'c');
% plot(aList, E0_mf, '--r');
% plot(aList, E0_exact, '--');
xlabel('a');
ylabel('Energy');
title('Ground State Energy Estimation');
% legend('CI - 2 orbs', 'CI - 5 orbs', 'CI - 10 orbs', 'MF', 'Exact');
legend('CI - 2 orbs', 'CI - 5 orbs', 'CI - 10 orbs');
hold off;
