%% Problem 1
clear all;
close all;
clc;

echodemo('h2o.m', 1);
echodemo('h2o.m', 2);
echodemo('h2o.m', 3);
echodemo('h2o.m', 4);

fprintf('H2O Ground State Energy: %.3f\n', energy);

%% Problem 2
clear all;
close all;
clc;

echodemo('h2o.m', 1);
echodemo('h2o.m', 2);
echodemo('h2o.m', 3);
echodemo('h2o.m', 4);

PS = P * S;
q_H1 = chg(1) - PS(1,1);
q_H2 = chg(2) - PS(2,2);
q_O = chg(3) - trace(PS(3:7, 3:7));

fprintf('Effective Charge (Mulliken Population Analysis):\n');
fprintf('H1 - %.3f\n', q_H1);
fprintf('H2 - %.3f\n', q_H2);
fprintf('O  - %.3f\n', q_O);

%% Problem 3
clear all;
close all;
clc;

len_list = 1.0:0.1:2.5;
len_energies = zeros(1, length(len_list));

for li = 1:length(len_list)
    echodemo('h2o.m', 1);
    len_OH1 = len_list(li);
    len_OH2 = len_list(li);
    echodemo('h2o.m', 2);
    echodemo('h2o.m', 3);
    echodemo('h2o.m', 4);
    len_energies(li) = energy;
end

plot(len_list, len_energies);
xlabel('Bond Lengths');
ylabel('Energy');

%% Problem 4
clear all;
close all;
clc;

ang_list = 80:1:105;
ang_energies = zeros(1, length(ang_list));

for li = 1:length(ang_list)
    echodemo('h2o.m', 1);
    ang_HOH = degtorad(ang_list(li));
    echodemo('h2o.m', 2);
    echodemo('h2o.m', 3);
    echodemo('h2o.m', 4);
    ang_energies(li) = energy;
end

plot(ang_list, ang_energies);
xlabel('Bond Angles (deg)');
ylabel('Energy');

%% Problem 5
clear all;
close all;
clc;

len_list = 1.995:0.001:2.005;
len_energies = zeros(1, length(len_list));

for li = 1:length(len_list)
    echodemo('h2o.m', 1);
    len_OH1 = len_list(li);
    echodemo('h2o.m', 2);
    echodemo('h2o.m', 3);
    echodemo('h2o.m', 4);
    len_energies(li) = energy;
end

plot(len_list, len_energies);
xlabel('Bond Lengths');
ylabel('Energy');

%% Problem 6
% clear all;
% close all;
% clc;

format long;
%echodemo('ps6.m', 5);
[sort_energs, index] = sort(len_energies);

delta = len_list(2) - len_list(1);

der2_E = len_energies(index(1) + 1) - 2*len_energies(index(1)) + len_energies(index(1) - 1);
der2_E = der2_E / delta^2

%% Problem 7

p_e_ratio = 1836.153;
m_OH = (1*p_e_ratio) * (9*p_e_ratio) / ((1*p_e_ratio) + (9*p_e_ratio));
w_OH = sqrt(der2_E / m_OH);

w_OH_wn = w_OH * 2.2 * 10^5

%% Problem 8
clear all;
close all;
clc;

echodemo('h2o.m', 1);
% ang_HOH = degtorad(95);
l = 1.85;
len_OH1 = l;
len_OH2 = l;    
echodemo('h2o.m', 2);
echodemo('h2o.m', 3);
echodemo('h2o.m', 4);

eigs