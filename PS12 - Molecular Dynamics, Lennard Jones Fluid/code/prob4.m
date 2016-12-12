%% Problem 4 part (v)
clear all;

N = 25;                 % number of particles
density = 0.9;    % density of particles in box
T = 0.45;                  % target temperature

dt=0.01;           % time interval per step
nsteps = 5000;     % number of steps to run simulation for

echodemo('LJSimulation', 2)

nmax = 150;
count = zeros(1, nmax+1);
corr = zeros(1, nmax+1);
for m = 1:nsteps-nmax-1000
    for n = 0:nmax
        count(n+1) = count(n+1) + 2 * length(vxtraj(m, :));
        corr(n+1) = corr(n+1) + sum(vxtraj(m, :) .* vxtraj(m + n, :));
        corr(n+1) = corr(n+1) + sum(vytraj(m, :) .* vytraj(m + n, :));
    end
end

corr = corr ./ count;
corr = corr / corr(1);

% plot(dt * (0:nmax), corr);
% title('\rho = 0.75');
% xlabel('t');
% ylabel('C(t)')
