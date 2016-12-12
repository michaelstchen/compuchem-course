%% Problem 2 part (ii)

N = 25;                 % number of particles
density = 0.5102041;    % density of particles in box
T = 1;                  % target temperature

dt=0.01;            % time interval per step
nsteps = 1000;      % number of steps to run simulation for

echodemo('LJSimulation.m', 2);

figure
hold on
plot(ttraj, Utraj)
plot(ttraj, Ktraj, 'r')
plot(ttraj, Etraj, 'g')
xlabel('time')
ylabel('Energy')

%% Problem 2 part (v)
clear all;
close all;
clc;

N = 25;                 % number of particles
density = 0.5102041;    % density of particles in box
T = 0.45;                  % target temperature

dt=0.01;            % time interval per step
nsteps = 5000;      % number of steps to run simulation for

echodemo('LJSimulation.m', 2);

figure
plot(ttraj, Tefftraj)
xlabel('time')
ylabel('T_{eff}')

%% Problem 2 part (vi)
Teff_avg = sum(Tefftraj(1001:end)) / length(Tefftraj(1001:end));
