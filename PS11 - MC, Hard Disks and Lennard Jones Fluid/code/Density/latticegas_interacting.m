clear

L=6;
M=L*L;

% mus = -3:0.1:-1;
% size = numel(mus);
% avden = zeros(1,size);
% for mui=1:size
%   mu = mus(mui)

mu=-2;
T=0.55;
beta=1/T;

n = zeros(L,L);
N=0;
E=0;

positions = zeros(2,M);

nsteps=100000;

avN = 0;

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
  if (mod(step,1000)==0)
    drawlatticeconfig(L,n);
    drawnow;
  end
end

% drawlatticeconfig(L,n);
avN = avN/nsteps;

% avden(mui) = avN/M;
% 
% end
% 
% clf
% plot(mus,avden)