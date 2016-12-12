% Problem 5
clear all;
close all;
clc;

K = 25;
n = (sqrt(K) - 1) / 2;
alpha = 2;
space = 0.5;

[energ, wf] = harmOsc2DEnergies(K, alpha, space);
[energ_sort, index] = sort(energ);

% % Ground state
% c = wf(:,index(0 + 1));
% One of the first excited state
c = wf(:, index(1 + 10));

xvals = -3:0.2:3;
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

surf(xvals, yvals, psi);
% title('Ground State 2D HO Wavefn');
title('1st Excited 2D HO Wavefn');
xlabel('x');
ylabel('y');
zlabel('Energy');
