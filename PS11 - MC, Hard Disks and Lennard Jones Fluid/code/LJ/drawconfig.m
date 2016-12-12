function drawconfig(positions,N,L,Lmax)

clf
hold on
for i=1:N
  r = positions(:,i);
  r = r -L*floor(r/L);
  plotcircle([r(1) r(2)],0.5);
end

x = [0 L L 0 0];
y = [0 0 L L 0];
plot(x,y,'k')

xlim([-0.5 Lmax]);
ylim([-0.5 Lmax]);

drawnow