%% Problem 1 part (ii)
clear all;
close all;
clc;

alpha=2;
alpha2=alpha^2;
% the spacing btw basis fn peaks
spacing=0.5;

% K is the total number of basis fns
K = 9;
n = (K - 1) / 2;

% S is out overlap matrix
S = zeros(K,K);
% H is our Hamiltonian matrix
H = zeros(K,K);

% Calculate matrix elements
for A=1:K
    xA = spacing*(-n+A-1);
    for B=1:K
        xB = spacing*(-n+B-1);
        
        deltax = xA-xB;
        deltax2 = deltax^2;
        
        sumx = xA+xB;
        sumx2 = sumx^2;
        
        S(A,B) = sqrt(pi/(2*alpha))*exp(-0.5*alpha*deltax2);
        H(A,B) = 0.5*S(A,B)*(alpha - alpha2*deltax2 ...
            + 0.25*(1/alpha + sumx2));
    end
end

% Solve for energies (eigenvalues) and
% wavefunction basis coeff (eigenvectors)
Sinv = inv(S);
SinvH = Sinv*H;
[U D] = eig(SinvH);
eigs = diag(D);
[eigsort index] = sort(eigs);

fprintf('First 9 HO Energy Levels\n');
for i = 1:9
    fprintf('%d : %.4f\n', i-1, eigsort(i));
end


%% Problem 1 part (iii)
clear all;
close all;
clc;

alpha=2;
alpha2=alpha^2;
% the spacing btw basis fn peaks
spacing=0.5;

KArray = 9:2:121;
energies = zeros(length(KArray), 9);
i = 1;
for K = KArray
    % K is the total number of basis fns
    n = (K - 1) / 2;

    % S is out overlap matrix
    S = zeros(K,K);
    % H is our Hamiltonian matrix
    H = zeros(K,K);

    % Calculate matrix elements
    for A=1:K
        xA = spacing*(-n+A-1);
        for B=1:K
            xB = spacing*(-n+B-1);

            deltax = xA-xB;
            deltax2 = deltax^2;

            sumx = xA+xB;
            sumx2 = sumx^2;

            S(A,B) = sqrt(pi/(2*alpha))*exp(-0.5*alpha*deltax2);
            H(A,B) = 0.5*S(A,B)*(alpha - alpha2*deltax2 ...
                + 0.25*(1/alpha + sumx2));
        end
    end

    % Solve for energies (eigenvalues) and
    % wavefunction basis coeff (eigenvectors)
    Sinv = inv(S);
    SinvH = Sinv*H;
    [U D] = eig(SinvH);
    eigs = diag(D);
    [eigsort index] = sort(eigs);
    
    energies(i, :) = transpose(eigsort(1:9));
    i = i+1;
end

plot(KArray, energies);
title('Estimated HO Energy Convergence');
xlabel('K');
ylabel('Energy (\epsilon)');


