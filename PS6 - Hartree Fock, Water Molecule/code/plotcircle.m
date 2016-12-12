function plotcircle(center, r)

twopi = 2*pi;
npts=100;

t=0:twopi/npts:twopi;
x = center(1) + r*cos(t);
y = center(2) + r*sin(t);

plot(x,y,'Color',[0 0 1]);