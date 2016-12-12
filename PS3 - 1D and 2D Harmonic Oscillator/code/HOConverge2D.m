%% Problem 4 part(i)
clear all;
close all;
clc;

alpha = 2;
space = 0.5;
K = 49;

[energ, wf] = harmOsc2DEnergies(K, alpha, space);
[energ_sort, indices] = sort(energ);

fprintf('First 45 2D Wavfn Energies\n');
for i = 1:45
    fprintf('%d & %.4f \\\\ \n', i, energ_sort(i));
end

%% Problem 4 part(ii)
clear all;
close all;
clc;

alpha = 2;
space = 0.5;

nArray = 3:9;
KArray = (2 * nArray + 1).^2;

energies = zeros(length(KArray), 45);
i = 1;
for K = KArray
    [energ, wf] = harmOsc2DEnergies(K, alpha, space);
    [energ_sort, indices] = sort(energ);
    energies(i, :) = transpose(energ_sort(1:45));
    i = i + 1;
end

plot(KArray, energies);
title('Estimated HO Energy Convergence');
xlabel('K');
ylabel('Energy (\epsilon)');

