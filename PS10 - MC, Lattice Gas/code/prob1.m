%% Problem 1 part(ii)
L = 10;
M = L * L;

mu_list = -3:0.1:3;
T=1;
beta=1/T;

avN_dev = @(B, mu, M) M * (exp(B*mu) / (1 + exp(B*mu))^2);

avN_dev_list = zeros(1, length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    avN_dev_list(i) = avN_dev(beta, mu, M);
end

plot(mu_list, avN_dev_list);
xlabel('\mu');
ylabel('<\deltaN^2>');
title('Exact');

%% Problem 1 part(iii)

mu_list = -3:0.2:3;
avN_dev_list1 = zeros(1, length(mu_list));
avN_dev_list2 = zeros(1, length(mu_list));
for i = 1:length(mu_list)
    echodemo('latticegas_noninteracting', 1);
    mu = mu_list(i);
    nsteps = 100000;
    echodemo('latticegas_noninteracting', 2);
    
    avN = avN/nsteps;
    avN_dev_list1(i) = sum((N_list - avN).^2) / nsteps;
    avN_dev_list2(i) = sum((N_list(1000:end) - avN).^2) / (nsteps-1000);
end

subplot(2,1,1)
plot(mu_list, avN_dev_list1);
xlabel('\mu');
ylabel('<\deltaN^2>');
title('Simulation');

subplot(2, 1, 2)
plot(mu_list, avN_dev_list2);
xlabel('\mu');
ylabel('<\deltaN^2>');
title('Simulation (Truncated)');