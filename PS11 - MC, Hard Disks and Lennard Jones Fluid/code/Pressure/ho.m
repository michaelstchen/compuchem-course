clear

x=1;
v=0;
F=-x;

dt=0.01;
T = 2*pi*3;
nsteps = round(T/dt);

xtraj = zeros(1,nsteps);
ttraj = dt*(1:nsteps);

for step=1:nsteps
  v = v + (dt/2)*F;
  x = x + dt*v;
  F = -x;
  v = v + (dt/2)*F;
  
  xtraj(step) = x;
  
end

clf
plot(ttraj,xtraj,'o')
exact = cos(ttraj);
hold on
plot(ttraj,exact)
f2f