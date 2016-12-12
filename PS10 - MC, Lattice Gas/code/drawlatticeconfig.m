function drawlatticeconfig(L,n)

clf
hold on

for x=1:L
  for y=1:L
    if (n(x,y)==1)
      plotcircle([x+1/2 y+1/2],0.4);
    end
  end
end

for x=1:L+1
  linex = [x x];
  liney = [1 L+1];
  plot(linex,liney,'k')
end

for y=1:L+1
  linex = [1 L+1];
  liney = [y y];
  plot(linex,liney,'k')
end

xlim([0 L+2]);
ylim([0 L+2]);

axis equal
