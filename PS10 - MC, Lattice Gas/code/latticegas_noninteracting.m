%% Setting up parameters
L=10;       % Number of lattice cells in 1D
M=L*L;      % Total number of latticce cells (2D)

mu=1;       % Chemical potentials (rev work to remove part from bath)
T=1;        % Temperature
beta=1/T;   % Beta (in units Kb)

n = zeros(L,L);     % Keeps track of whether cell has particle
N=0;                % Current number of particles
avN = 0;            % Average number of particles 

positions = zeros(2,M);     % Positions of all particles

nsteps=10000;               % Number of simulation steps

%% Run MC Simulation

N_list = zeros(1, nsteps);  % Keeps track of N at each step
for step=1:nsteps
  if (rand < 0.5)
    x = ceil(rand*L);
    y = ceil(rand*L);
    if (n(x,y)==0 && rand < exp(beta*mu)*M/(N+1))
      % Add a particle
      n(x,y)=1;
      N = N+1;
      positions(1,N)=x;
      positions(2,N)=y;
    end
  elseif (N>0)
    if (rand < exp(-beta*mu)*N/M)
      % Remove a particle
      index = ceil(rand*N);
      x = positions(1,index);
      y = positions(2,index);
      n(x,y)=0;
      positions(1,index) = positions(1,N);
      positions(2,index) = positions(2,N);
      N=N-1;
    end
  end
  avN = avN + N;
  N_list(step) = N;
end

%% Drawing Lattice
drawlatticeconfig(L,n);
avN = avN/nsteps;

avN/M

z=exp(beta*mu);
z/(1+z)