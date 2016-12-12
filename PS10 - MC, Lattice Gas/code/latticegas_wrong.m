clear

L=10;
M=L*L;

n = zeros(L,L);
N=0;

positions = zeros(2,M);

nsteps=1000;

avN = 0;

for step=1:nsteps
  if (rand < 0.5)
    x = ceil(rand*L);
    y = ceil(rand*L);
    if (n(x,y)==0)
      n(x,y)=1;
      N = N+1;
      positions(1,N)=x;
      positions(2,N)=y;
    end
  elseif (N>0)
    if (rand < N/M)
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
end

drawlatticeconfig(L,n);
avN = avN/nsteps;

avN