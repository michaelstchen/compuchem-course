%% Problem 2 part (ii)
clear all;
close all;
clc;

rcut = 2.5;

rlist = 0.95:0.01:3.0;
ulist = zeros(1, length(rlist));

for iu = 1:length(ulist)
    ulist(iu) = LJPotential(rlist(iu), rcut);
end

plot(rlist, ulist)
xlabel('r')
ylabel('u(r)')

%% Problem 2 part (iv)
clear all;
close all;
clc;

Nside=5;

N = Nside^2;
Llist = 5:0.1:10;
Ulist = zeros(1, length(Llist));
r = zeros(2,N);

for li = 1:length(Llist)
    L = Llist(li);
    spacing = L/Nside;
    index = 0;
    for i=1:Nside
      for j=1:Nside
        index = index + 1;
        x = i*spacing+spacing/2;
        y = j*spacing+spacing/2;
        r(:,index) = [x y]';
      end
    end

    Ulist(li) = LJPotentialTotal(N, L, r);
end

plot(Llist, Ulist)
xlabel('L')
ylabel('U')

%% Problem 2 part (v)
clear all;
close all;
clc;

Nside=5;

N = Nside^2;
L = 7;
r = zeros(2,N);

spacing = L/Nside;
index = 0;
for i=1:Nside
  for j=1:Nside
    index = index + 1;
    x = i*spacing+spacing/2;
    y = j*spacing+spacing/2;
    r(:,index) = [x y]';
  end
end

Uorig = LJPotentialTotal(N, L, r);

r(1, 1) = r(1, 1) + 0.05;
r(2, 1) = r(2, 1) + 0.0866;

Ushift = LJPotentialTotal(N, L, r);

fprintf('Change in Energy = %.3f\n', Ushift - Uorig);

%% Problem 2 part (vi) (vii)
clear all;
close all;
clc;

echodemo('LJSimulation.m', 1)
echodemo('LJSimulation.m', 2)

figure
plot(hr.vals,hr.hist./(nsweeps*pi*N*density*hr.binwidth*hr.vals))
xlabel('r')
ylabel('g(r)')

figure
plot(hu.vals, hu.hist/hu.count)
xlabel('U')
ylabel('p(U)')
