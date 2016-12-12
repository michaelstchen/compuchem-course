%% Parameter setup
L=10;
M=L*L;

mu=-2;
T=0.5;
beta=1/T;

nsteps=100000;

h.count=0;
h.range=[-0.5 M+0.5];
h.binwidth=1; 

%% Running MC Simulation
n = zeros(L,L);
avN = 0;
N=0;
E=0;

positions = zeros(2,M);

for step=1:nsteps
  if (rand < 0.5)
    x = ceil(rand*L);
    y = ceil(rand*L);
    nn = nneigh(L,n,x,y);
    deltaE = -nn;
    if (n(x,y)==0 && rand < exp(beta*mu-beta*deltaE)*M/(N+1))
      n(x,y)=1;
      N = N+1;
      E = E+deltaE;
      positions(1,N)=x;
      positions(2,N)=y;
    end
  elseif (N>0)
    index = ceil(rand*N);
    x = positions(1,index);
    y = positions(2,index);
    nn = nneigh(L,n,x,y);
    deltaE = nn;
    if (rand < exp(-beta*mu-beta*deltaE)*N/M)
      n(x,y)=0;
      positions(1,index) = positions(1,N);
      positions(2,index) = positions(2,N);
      N=N-1;
      E = E+deltaE;
    end
  end
  avN = avN + N;
  h = histo(h, N);

end

%% Draw Lattice Config
drawlatticeconfig(L,n);
avN = avN/nsteps;

avN/M
