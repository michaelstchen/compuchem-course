%% Initializations

T = 1;          % dimensionless temperature
d = 1;          % magnitude/weight of displacement

x = 0;          % current position
u = 0.5*x^2;    % current energy

% histogram initializations
h.count=0;
h.range=[-3 3];
h.binwidth=0.1;

nsteps=100000;              % number of iterations to run MC simulation
num_accepted = 0;           % number of accepted moves

%% MC simulation for-loop
xtraj = zeros(1,nsteps);    % list of positions at each iteration
for step=1:nsteps
  xtrial = x + (rand-0.5)*d;        % trial displacement
  utrial = 0.5*xtrial^2;            % trial energy
  
  % Check if we should accept move (Glauber)
  if ( rand < 1 / ( 1 + exp((utrial-u)/T) ) )
    x = xtrial;
    u = utrial;
    num_accepted = num_accepted + 1;
  end
  
  % update statistics
  xtraj(step)=x;
  h = histo(h,x);
end

%% Plotting results

hold on
% plot(1:nsteps,xtraj)
plot(h.vals,h.hist/(h.count*h.binwidth))
plot(h.vals,exp(-0.5*(h.vals .^2)/T)/sqrt(2*pi*T))