%% Problem 1 part (ii)
clear all;
close all;
clc;

alpha=2;
alpha2=alpha^2;     
spacing=0.5;

n=5;
K=2*n+1;

S = zeros(K,K);
h = zeros(K,K);
G = zeros(K,K);

for A=1:K
    xA = spacing*(-n+A-1);
    for B=1:K
        xB = spacing*(-n+B-1);
        
        S(A,B) = s(xA,xB,alpha);
        h(A,B) = f(xA,xB,alpha);
        G(A,B) = g(xA,xB,alpha);
        
    end
end


Sinv = inv(S);
Sinvh = Sinv*h;
[U D] = eig(Sinvh);
eigs = diag(D);
[eigs index] = sort(eigs);

c = U(:,index(1));
norm = 1/sqrt(c'*S*c);
c = c*norm;

aList = 0:0.00025:0.01;
mf_evals = zeros(1, length(aList));
pert_evals = zeros(1, length(aList));
niter=10;

for a = 1:length(aList)
    for iter=1:niter
        avy4 = c'*G*c;
        heff = h + aList(a)*avy4*G;

        Sinvheff = Sinv*heff;
        [U D] = eig(Sinvheff);
        eigs = diag(D);
        [eigs index] = sort(eigs);

        c = U(:,index(1));
        norm = 1/sqrt(c'*S*c);
        c = c*norm;

        energy = 2*c'*h*c + aList(a)*(c'*G*c)^2;

    end
    mf_evals(a) = energy;
    pert_evals(a) = 1 + 9*aList(a)/16;
end


hold on;
plot(aList, mf_evals);
plot(aList, pert_evals, 'o');
xlabel('a')
ylabel('Energy')
title('Ground State Energies from Mean Field and Perturbation Theories')
legend('mean field', 'pertubation');

%% Problem 1 part (iii) and (iv)
clear all;
close all;
clc;

alpha=2;
alpha2=alpha^2;     
spacing=0.5;

n=5;
K=2*n+1;

S = zeros(K,K);
h = zeros(K,K);
G = zeros(K,K);

for A=1:K
    xA = spacing*(-n+A-1);
    for B=1:K
        xB = spacing*(-n+B-1);
        
        S(A,B) = s(xA,xB,alpha);
        h(A,B) = f(xA,xB,alpha);
        G(A,B) = g(xA,xB,alpha);
        
    end
end


Sinv = inv(S);
Sinvh = Sinv*h;
[U D] = eig(Sinvh);
eigs = diag(D);
[eigs index] = sort(eigs);

c = U(:,index(1));
norm = 1/sqrt(c'*S*c);
c = c*norm;

aList = 0:100;
mf_evals = zeros(1, length(aList));
niter=10;

for a = 1:length(aList)
    for iter=1:niter
        avy4 = c'*G*c;
        heff = h + aList(a)*avy4*G;

        Sinvheff = Sinv*heff;
        [U D] = eig(Sinvheff);
        eigs = diag(D);
        [eigs index] = sort(eigs);

        c = U(:,index(1));
        norm = 1/sqrt(c'*S*c);
        c = c*norm;

        energy = 2*c'*h*c + aList(a)*(c'*G*c)^2;

    end
    mf_evals(a) = energy;
end

K = 225;
E0 = zeros(1, 100);
for a = 1:length(aList)
    [energ, wf] = HONonLin2D(K, alpha, spacing, aList(a));
    [energ_sort, indices] = sort(energ);
    E0(a) = energ_sort(1);
end

figure;
plot(aList, mf_evals);
xlabel('a')
ylabel('Energy')
title('Ground State Energies from Mean Field Theory')
legend('mean field', 'exact')

figure;
hold on;
plot(aList, mf_evals);
plot(aList, E0, 'o');
xlabel('a')
ylabel('Energy')
title('Comparing MFT with Exact Values')
legend('mean field', 'exact')
hold off;


%% Problem 1 part (v)
clear all;
close all;
clc;

alpha=2;
alpha2=alpha^2;     
spacing=0.5;

n=5;
K=2*n+1;

S = zeros(K,K);
h = zeros(K,K);
G = zeros(K,K);

for A=1:K
    xA = spacing*(-n+A-1);
    for B=1:K
        xB = spacing*(-n+B-1);
        
        S(A,B) = s(xA,xB,alpha);
        h(A,B) = f(xA,xB,alpha);
        G(A,B) = g(xA,xB,alpha);
    end
end


Sinv = inv(S);
Sinvh = Sinv*h;
[U D] = eig(Sinvh);
eigs = diag(D);
[eigs index] = sort(eigs);

c = U(:,index(1));
norm = 1/sqrt(c'*S*c);
c = c*norm;

a = 10;
niter=10;

for iter=1:niter
    avy4 = c'*G*c;
    heff = h + a*avy4*G;

    Sinvheff = Sinv*heff;
    [U D] = eig(Sinvheff);
    eigs = diag(D);
    [eigs index] = sort(eigs);

    c = U(:,index(1));
    norm = 1/sqrt(c'*S*c);
    c = c*norm;
end

xvals = -2:0.1:2;
yvals = xvals;
psi = zeros(length(xvals), length(yvals));
for ix = 1:length(xvals)
    pA = 0;
    for i = 1:K
        xA = spacing*(-n+i-1);
        pA = pA + c(i) * exp(-alpha * (xvals(ix) - xA)^2);
    end

    for iy = 1:length(yvals)
        pB = 0;
        for i = 1:K
            xB = spacing*(-n+i-1);
            pB = pB + c(i) * exp(-alpha * (yvals(iy) - xB)^2);
        end

        psi(ix, iy) =  pA * pB;
    end
end

i_mid = round(length(psi)/2);
psi = psi * sign(psi(i_mid, i_mid)); 

contourm(xvals, yvals, psi);
title('Ground State Contour Maps');
xlabel('x');
xlim([-2 2]);
ylabel('y');
ylim([-2 2]);
axis equal