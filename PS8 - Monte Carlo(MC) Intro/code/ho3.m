% Used for Problem 2 parts (ix) and (x)
clear all;
close all;

M=30;
E_list = 1:30;
beta = zeros(1, length(E_list));
for iE = 1:length(E_list)
    E = E_list(iE); 
    n = zeros(1,M);

    n(5)=E;

    nsteps = 10000;

    h.count=0;
    h.range = [-0.5, E+0.5];
    h.binwidth = 1;
    for step=1:nsteps
      i = ceil(rand*M);
      j = ceil(rand*M);

      if (n(i)>0)
        n(i) = n(i)-1;
        n(j) = n(j)+1;
      end

      for i = 1:length(n)
          h = histo(h,n(i));
      end
    end
    
    prob = h.hist / h.count;
    beta(iE) = -log(prob(2) / prob(1));
end

prob = h.hist / h.count;

% % FOR PROBLEM 2 PART 9
% beta = -log(prob(2) / prob(1));
% alpha = log(prob(1));
% x = [0, M/2];
% y = [alpha, alpha-beta*M/2];
% hold on;
% plot(0:M, log(prob));
% plot(x, y, 'r');
% title('Exponential Probability Decay')
% xlabel('n');
% ylabel('ln[p(n)]');

% FOR PROBLEM 2 PART 10
scatter(1 ./ beta, E_list / M);
title('Beta for Simulations with Various E');
ylabel('Expected n');
xlabel('1/Beta');

fprintf('beta & simu n & known n & error \\\\ \n');
for i = 1:length(beta)
    beta_i = beta(i);
    sim_ni = E_list(i) / M;
    known_ni = 1 / (exp(beta_i) - 1);
    per_error = abs(known_ni-sim_ni) / known_ni;
    fprintf('%.3f & %.3f & %.3f & %.3f \\\\ \n', beta_i, sim_ni, known_ni, per_error);
end