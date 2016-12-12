% Used for Problem 2 parts (ii) and (iii)
clear all;
close all;

M=100;               % Number of oscillators
E=500;               % Total amount of energy

n = zeros(1,M);     % Distribution of microstates

% Arbitarily initialize distrib
for i = 1:length(n)
    n(i) = E / M;
end

nsteps = 10000;    % Number of simulation iterations

% initialize histogram
h.count=0;
h.range = [-0.5, E+0.5];
h.binwidth = 1;

% stores data for part (ii) and (iii)
n1_vals = zeros(1, nsteps);
n1_means = [];

% Simulation
for step=1:nsteps
  % Randomly select two of our oscillators
  i = ceil(rand*M);
  j = ceil(rand*M);
  
  % Move a quanta from osc I (if quanta > 0) to osc J
  if (n(i)>0)
    n(i) = n(i)-1;
    n(j) = n(j)+1;
  end
   
  % update histogram
  h = histo(h,n(1));
  n1_vals(step) = n(1);
end

figure;
% plot normalized histogram
bar(h.vals,h.hist/h.count)

avg = 0;
for i = 1:length(h.hist)
    avg = avg + (h.hist(i)/h.count) * h.vals(i);
end

figure;
% plots for parts (ii) and (iii)
plot(1:nsteps, n1_vals)
title('Number of quanta in oscillator #1');
xlabel('Step')
ylabel('Count')