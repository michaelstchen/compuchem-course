clear all;
close all;
clc;

%%Q2 (ii)
x1 = 0:0.01:1;
y1 = [];
figure;
subplot(2,1,1);
hold on;
for i = 1:4
    for j = 0:0.01:1
        y1(end+1) = sin(pi * i * j);
    end
    plot(x1, y1);
    y1 = [];
end
hold off;

%% Q4 (iv)
a = 8 / (9 * pi);
gaussVar = @(x) (pi / (2 * a))^-0.75 * exp(-a * x^2);
realWF = @(x) pi^-0.5 * exp(-x);
x2 = 0:0.1:4;
y2_gauss = [];
y2_real = [];
for i = 0:0.1:4
    y2_gauss(end+1) = gaussVar(i);
    y2_real(end+1) = realWF(i);
end

subplot(2,1,2);
hold on;
plot(x2, y2_gauss, '-');
plot(x2, y2_real);
ylim([0 0.5]);
hold off;