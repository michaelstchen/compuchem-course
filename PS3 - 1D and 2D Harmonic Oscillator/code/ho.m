clear all;
close all;
clc;

% The width parameter of our Gaussian basis fns
alpha=2;
% The squared value of alpha (saves computation time)
alpha2=alpha^2;
% The spacing between the peaks of our basis fns
spacing=0.5;

% K is total number of basis functions, calculated from
% an integer value n
n=20;
K=2*n+1;

% S is our matrix of overlap matrix elements
S = zeros(K,K);
% H is our matrix of Hamiltonian matrix elements
H = zeros(K,K);

% Calculating the matrix elements for S and H
% for every pair permutation of basis fns
for A=1:K
    % centering of basis fn A for this iteration
    xA = spacing*(-n+A-1);
    for B=1:K
        % centering of basis fn B for this iteration
        xB = spacing*(-n+B-1);
        
        % seperation btw the two bases for this iteration
        deltax = xA-xB;
        % the square of the separation
        deltax2 = deltax^2;
        
        % sum of the centered positions (for avg calc)
        sumx = xA+xB;
        % square of the sum
        sumx2 = sumx^2;
        
        S(A,B) = sqrt(pi/(2*alpha))*exp(-0.5*alpha*deltax2);
        H(A,B) = 0.5*S(A,B)*(alpha - alpha2*deltax2 ...
            + 0.25*(1/alpha + sumx2));
    end
end

% Finding energies (eigenvalues) and basis coeff (eigenvectors)
Sinv = inv(S);                  % inverse of the overlap matrix S
SinvH = Sinv*H;                 % the value we get from left multiplying H by the inverse of S
[U D] = eig(SinvH);             % the columns of U are our basis coeff, the diagonals of D are our energies
eigs = diag(D);                 % a list representation of the diagonals/energies from D
[eigsort index] = sort(eigs);   % a sorted list of energies with a list (index) mapping the current and former indices

% c is a column vector of basis coeffs for the nth sorted wavefunction
% based on energy where n is set in 'index(n)'
c = U(:,index(10));
% norm is the normalization constant for the wavefn specified by the
% basis coeffs in c
norm = 1/sqrt(c'*S*c);
c = c*norm;

% xvals is a list of x values we will plot for our HO wavefunction
xvals = -6:0.01:6;
% psi is a list of the y vals for our estimated wavefn
psi = 0*xvals;

% constructing our estimated wavefunction (summing up all of the
% component basis fns with their appropriate weights)
for A=1:K
    % xA is the centering for the basis fn examined in this iteration
    xA = spacing*(-n+A-1);
    psi = psi + c(A)*exp(-alpha*(xvals-xA).^2);
end
psi = psi*sign(psi(1));

% exact theoretical ground state yvals
psiexact = (pi^(-0.25))*exp(- 0.5*xvals .^ 2);
% exact theoretical n = 9 yvals
psiexact = psiexact .* (512*xvals.^9-9216*xvals.^7 ...
    +48384*xvals.^5 -80640*xvals.^3 + ...
    30240*xvals)/(2304*sqrt(35));
psiexact = psiexact*sign(psiexact(1));

% plotting our results
hold on
plot(xvals,psi,'r')
plot(xvals,psiexact, 'b')







