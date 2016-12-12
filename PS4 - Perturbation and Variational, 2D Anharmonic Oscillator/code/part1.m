%% Problem 1 part(ii)
clear all;
close all;
clc;

alpha = 2;
space = 0.5;
K = 225;

fprintf('Variational Ground State Energies\n');
for a = [0.001, 0.01, 0.1, 1, 10, 100]
    [energ, wf] = HONonLin2D(K, alpha, space, a);
    [energ_sort, indices] = sort(energ);
    fprintf('%.3f & %.4f \\\\ \n', a, energ_sort(1));
end
fprintf('\n');

%% Problem 1 part (iii)
clear all;
close all;
clc;

alpha = 2;
E0 = 1;
E_anharm = 9 / (16 * alpha^4);

fprintf('Perturbation Ground State Energies\n');
for a = [0.001, 0.01, 0.1, 1, 10, 100]
    fprintf('%.3f & %.4f \\\\ \n', a, E0 + a * E_anharm);
end

%% Problem 1 part (iv)
clear all;
close all;
clc;

alpha = 2;
space = 0.5;
K = 49;

E0 = zeros(1, 100);
aList = 1:100;
for a = aList
    [energ, wf] = HONonLin2D(K, alpha, space, a);
    [energ_sort, indices] = sort(energ);
    E0(a) = energ_sort(1);
end

plot(aList, E0);
title('Variational Ground State Energy Estimates for Various Anharmonicity Factors');
xlabel('Anharmonicity Factor (a)');
ylabel('Energy');


%% Problem 1 part (v) and (vi)
clear all;
close all;
clc;

alpha = 2;
space = 0.5;
K = 225;
n = (sqrt(K) - 1) / 2;

aList = [0.001, 0.01, 0.1, 1, 10, 100];
for ai = 1:length(aList)
    a = aList(ai);
    [energ, wf] = HONonLin2D(K, alpha, space, a);
    [energ_sort, indices] = sort(energ);
    
    c = wf(:, indices(1));

    xvals = -2:0.1:2;
    yvals = xvals;
    psi = zeros(length(xvals), length(yvals));
    for ix = 1:length(xvals)
        for iy = 1:length(yvals)
            i = 1;
            for xA = -n*space:space:n*space
                for yA = -n*space:space:n*space
                    psi(ix, iy) = psi(ix, iy)...
                                + c(i)*exp(-alpha*(xvals(ix) - xA).^2 ...
                                           - alpha*(yvals(iy) - yA).^2);
                    i = i + 1;
                 end
            end
        end
    end

    i_mid = round(length(psi)/2);
    psi = psi * sign(psi(i_mid, i_mid));    
    subplot(2, 3, ai);
    contourm(xvals, yvals, psi);
    title(sprintf('Ground State WaveFunc for A=%.3f', a));
    xlabel('x');
    xlim([-2 2]);
    ylabel('y');
    ylim([-2 2]);
    zlabel('Energy');
end