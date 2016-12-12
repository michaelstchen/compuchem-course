%% Problem 3 part (ii)
% clear all;
% close all;
% clc;

alpha=2;
alpha2=alpha^2;
spacing=0.5;

n=5;
K=2*n+1;

N=2;

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

c=zeros(K,N);
for j=1:N
    c(:,j) = U(:,index(j));
    norm = 1/sqrt(c(:,j)'*S*c(:,j));
    c(:,j) = c(:,j)*norm;
end
P = c*c';

aList = 0:0.2:10;
E0_mf = zeros(1, length(aList));
niter=10;

for ia = 1:length(aList)
    a = aList(ia);
    for iter=1:niter
        heff = h + a*( (trace(P*G))*G - G*P*G );

        Sinvheff = Sinv*heff;
        [U D] = eig(Sinvheff);
        eigs = diag(D);
        [eigs index] = sort(eigs);

        for j=1:N
            c(:,j) = U(:,index(j));
            norm = 1/sqrt(c(:,j)'*S*c(:,j));
            c(:,j) = c(:,j)*norm;
        end
        P = c*c';
        energy = trace(P*h) +0.5*a*((trace(P*G))^2 - trace(P*G*P*G) );
    end
    E0_mf(ia) = energy;
end

echodemo('prob2.m', 1);
E0_exact = zeros(1, length(aList));
for i = 1:length(aList)
    a = aList(i);
    [wf, D] = eig(S \ (Hx + Hy + a*A));
    energ = diag(D);
    [energ_sort, indices] = sort(energ);
    E0_exact(i) = energ_sort(1);
end

figure;
hold on;
plot(aList, E0_mf);
plot(aList, E0_exact, 'o');
xlabel('a');
ylabel('Energy');
title('Ground State Energy Estimation');
legend('mean field', 'exact');

% %% Problem 3 part (iii)
% clear all;
% close all;
% clc;
% 
% alpha=2;
% alpha2=alpha^2;
% spacing=0.5;
% 
% n=5;
% K=2*n+1;
% 
% N=2;
% 
% S = zeros(K,K);
% h = zeros(K,K);
% G = zeros(K,K);
% 
% for A=1:K
%     xA = spacing*(-n+A-1);
%     for B=1:K
%         xB = spacing*(-n+B-1);
%         
%         S(A,B) = s(xA,xB,alpha);
%         h(A,B) = f(xA,xB,alpha);
%         G(A,B) = g(xA,xB,alpha);
%         
%     end
% end
% 
% Sinv = inv(S);
% Sinvh = Sinv*h;
% [U D] = eig(Sinvh);
% eigs = diag(D);
% [eigs index] = sort(eigs);
% 
% c=zeros(K,N);
% for j=1:N
%     c(:,j) = U(:,index(j));
%     norm = 1/sqrt(c(:,j)'*S*c(:,j));
%     c(:,j) = c(:,j)*norm;
% end
% P = c*c';
% 
% a = 10;
% niter=10;
% 
% for iter=1:niter
%     heff = h + a*( (trace(P*G))*G - G*P*G );
% 
%     Sinvheff = Sinv*heff;
%     [U D] = eig(Sinvheff);
%     eigs = diag(D);
%     [eigs index] = sort(eigs);
% 
%     for j=1:N
%         c(:,j) = U(:,index(j));
%         norm = 1/sqrt(c(:,j)'*S*c(:,j));
%         c(:,j) = c(:,j)*norm;
%     end
%     P = c*c';
%     energy = trace(P*h) +0.5*a*((trace(P*G))^2 - trace(P*G*P*G) );
% end
% 
% c1 = U(:, index(1));
% c2 = U(:, index(2));
% 
% xvals = -2:0.1:2;
% yvals = xvals;
% psi = zeros(length(xvals), length(yvals));
% for ix = 1:length(xvals)
%     for iy = 1:length(yvals)
%         X1_x = 0;
%         X2_x = 0;
%         X1_y = 0;
%         X2_y = 0;
%         for i = 1:K
%             xA = spacing*(-n+i-1);
%             X1_x = X1_x + c1(i) * exp(-alpha * (xvals(ix) - xA)^2);
%             X2_x = X2_x + c2(i) * exp(-alpha * (xvals(ix) - xA)^2);
%             
%             X1_y = X1_y + c1(i) * exp(-alpha * (yvals(iy) - xA)^2);
%             X2_y = X2_y + c2(i) * exp(-alpha * (yvals(iy) - xA)^2);
%         end
% 
%         psi(ix, iy) = (X1_x*X2_y - X2_x*X1_y) / sqrt(2);
%     end
% end
% 
% contourm(xvals, yvals, psi, 20);
% xlabel('x');
% ylabel('y');
% title('Mean Field Antisymmetric Oscillator Ground State Wavefn');