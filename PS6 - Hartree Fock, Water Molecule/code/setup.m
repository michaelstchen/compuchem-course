clear all;
close all;
clc;

N=2;
Nnuc=2;
K=2;
L=3;

zeta = [2.0925 1.24];

alpha = zeros(L,K);
d = zeros(L,K);

alpha(:,1) = [0.109818 0.405771 2.22766]*zeta(1)^2;
alpha(:,2) = [0.109818 0.405771 2.22766]*zeta(2)^2;

d(:,1) = [0.444635 0.535328 0.154329]*(2/pi)^(3/4);
d(:,2) = [0.444635 0.535328 0.154329]*(2/pi)^(3/4);

d = d .* alpha.^(3/4);

dist = 1.4632;
Rnuc = zeros(3,Nnuc);
Rnuc(1,2) = dist;
R = zeros(3,L,K);
for p=1:L
  R(1,p,2) = dist;
end

xvals = -1:0.1:3;
yvals = -2:0.1:2;
sz = numel(xvals);
psi = zeros(sz,sz);

for i=1:sz
  x = xvals(i);
  for j=1:sz
    y = yvals(j);
    for p=1:L
      dx = x - R(1,p,1);
      dy = y - R(2,p,1);
      psi(j,i) = psi(j,i) + d(p,1)*exp(-alpha(p,1)*(dx^2+dy^2));
    end
  end
end

contour(xvals,yvals,psi,10)
axis equal
hold on

plotcircle([0 0 0],0.2)
plotcircle([dist 0 0],0.1)

% clf
% rvals = 0:0.01:4;
% psiex = sqrt(zeta(1)^3/pi)*exp(-rvals*zeta(1));
% plot(rvals,psiex)
% psi = 0*rvals;
% for p=1:L
%   psi = psi + d(p,1)*exp(-alpha(p,1)*rvals.^2)
% end
% hold on
% plot(rvals,psi)



