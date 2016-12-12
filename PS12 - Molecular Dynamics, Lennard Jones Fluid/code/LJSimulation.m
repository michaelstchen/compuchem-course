%% SETUP PARAMETERS
N = 25;                 % number of particles
density = 0.75;    % density of particles in box
T = 1.2;                  % target temperature

dt=0.01;           % time interval per step
nsteps = 5000;     % number of steps to run simulation for

%% SIMULATION
L = sqrt(N/density);    % box dimension
r = zeros(2,N);         % particle positions
v = zeros(2,N);         % particle velociies

% Initialize particle configurations and velocities
spacing = L/sqrt(N);
index=0;
for i=1:sqrt(N)
  for j=1:sqrt(N)
    % Positions
    x = i*spacing + spacing/2;
    y = j*spacing + spacing/2;
    
    % Randomize velocities
    vx = rand - 1/2;
    vy = rand - 1/2;
    
    index = index+1;
    r(:,index) = [x y]';
    v(:,index) = [vx vy]';
  end
end

% Subtracting center-of-mass velocity
vavg = sum(v, 2) / N;
for i = 1:N
    v(:, i) = v(:, i) - vavg;
end

% Drawing initial configuration
% drawconfig(r,N,L,L+0.5);

% initialize histogram keeping track of distances
hr.count=0;
hr.range=[0, L/2];
hr.binwidth=0.1;
hr.hist = [];

% Initial Forces and Energies
[F, U, hr] = LJForcePot(r, 2.5, L, hr);
v2 = v.^2;
K=0.5*sum(v2(:));

Utraj = zeros(1,nsteps);    % keeps track of pot energy over simulation
Ktraj = zeros(1,nsteps);    % keeps track of kin energy over simulation
ttraj = dt*(1:nsteps);      % keeps track of positions over simulation
Tefftraj = zeros(1,nsteps); % keeps track of effective temperature
vxtraj = zeros(nsteps-1000,N);
vytraj = zeros(nsteps-1000,N);

% MD Simulation Loop
for step=1:nsteps
    % Verlat updates to velocities and positions
    v = v + (dt/2)*F;
    r = r + dt*v;
    [F, U, hr] = LJForcePot(r, 2.5, L, hr);
    v = v + (dt/2)*F;
    
    v2 = v.^2;
    K=0.5*sum(v2(:));

    Utraj(step) = U;
    Ktraj(step) = K;

    % Rescaling velocities to achieve T
    Teff = K / (N - 1);
    Tefftraj(step) = Teff;
    if (step <= 1000 && mod(step, 100) == 0)
        alpha = sqrt(T / Teff);
        v = alpha * v;
    end
    
    if (step > 1000)
        vxtraj(step-1000, :) = v(1, :);
        vytraj(step-1000, :) = v(2, :);
    end
%     if (mod(step,10)==0)
%         drawconfig(r,N,L,L+0.5);
%     end
  
end