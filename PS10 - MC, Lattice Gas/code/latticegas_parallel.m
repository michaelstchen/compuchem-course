%% Setup parameters
L = 10;
M = L*L;

mu = -2;
T = 0.4:0.1:1.0;
beta=1./T;
nreps = length(T);

n = zeros(L,L,nreps);
N=zeros(1,nreps);
E=zeros(1,nreps);
positions = zeros(2,M,nreps);

nsteps=10^6;

hist_list = struct();
for i = 1:nreps
    h_list(i).count = 0;
    h_list(i).range=[-0.5 M+0.5];
    h_list(i).binwidth=1; 
    h_list(i).numbins=101;
    h_list(i).hist = [];
    h_list(i).vals = [];
end

%% Run MC Simulation
for step=1:nsteps
  
  rep1 = ceil(rand*nreps);
  rep2 = ceil(rand*nreps);
  
  dN = N(rep2) - N(rep1);
  dE = E(rep2) - E(rep1);
  dE_tot = mu * dN - dE;
  dB = beta(rep2) - beta(rep1);
  
  if (rand < exp(dB * dE_tot))
    [E(rep1), E(rep2)] = deal(E(rep2),E(rep1));
    [N(rep1), N(rep2)] = deal(N(rep2),N(rep1));
    [n(:,:,rep1), n(:,:,rep2)] = deal(n(:,:,rep2),n(:,:,rep1));
    [positions(:,:,rep1), positions(:,:,rep2)] = ...
    deal(positions(:,:,rep2),positions(:,:,rep1));
  end
    
  for rep=1:nreps
    if (rand < 0.5)
      x = ceil(rand*L);
      y = ceil(rand*L);
      nn = nneigh(L,n(:,:,rep),x,y);
      deltaE = -nn;
      if (n(x,y,rep)==0 && ...
          rand < exp(beta(rep)*mu-beta(rep)*deltaE)*M/(N(rep)+1))
        n(x,y,rep)=1;
        N(rep) = N(rep)+1;
        E(rep) = E(rep)+deltaE;
        positions(1,N(rep),rep)=x;
        positions(2,N(rep),rep)=y;
      end
    elseif (N(rep)>0)
      index = ceil(rand*N(rep));
      x = positions(1,index,rep);
      y = positions(2,index,rep);
      nn = nneigh(L,n(:,:,rep),x,y);
      deltaE = nn;
      if (rand < exp(-beta(rep)*mu-beta(rep)*deltaE)*N(rep)/M)
        n(x,y,rep)=0;
        positions(1,index,rep) = positions(1,N(rep),rep);
        positions(2,index,rep) = positions(2,N(rep),rep);
        N(rep)=N(rep)-1;
        E(rep) = E(rep)+deltaE;
      end
    end
  end

  for i = 1:nreps
    h_list(i) = histo(h_list(i), N(i));
  end
end

%% Plotting results
for i = 1:length(T)
    subplot(3, 3, i)
    h = h_list(i);
    plot(h.vals, h.hist/h.count);
    title(sprintf('T = %.1f', T(i)));
    xlabel('N');
    xlim([0 M])
    ylabel('P(N)');
end