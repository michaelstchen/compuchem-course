clear
clf

% alpha sets the basis peaks' widths
alpha=2;
alpha2=alpha^2;     
% spacing between basis peaks
spacing=0.5;

n=5;
% K is the number of basis functions
K=2*n+1;

% S is the overlap matrix
S = zeros(K,K);
% H is our harmonic-oscillator Hamiltonian matrix
h = zeros(K,K);
% G is the matrix that accounts for the anharmonic elements
G = zeros(K,K);

% Calculating the matrix elements for S, h, G for every
% permutation of the K basis function pairs
for A=1:K
    % first peak center
    xA = spacing*(-n+A-1);
    for B=1:K
        % second peak center
        xB = spacing*(-n+B-1);
        
        S(A,B) = s(xA,xB,alpha);    % overlap matrix element
        h(A,B) = f(xA,xB,alpha);    % harmonic Hamiltonian matrix element
        G(A,B) = g(xA,xB,alpha);    % anharmonic Hamiltonian matrix element
        
    end
end

% Finding energies (eigenvalues) and basis coeff (eigenvectors)
% for our initial guess
Sinv = inv(S);
Sinvh = Sinv*h;
[U D] = eig(Sinvh);
eigs = diag(D);
[eigs index] = sort(eigs);
eigscopy= eigs;

% c is a normalized column vector of basis coeffs for the 
% ground state wavefunction
c = U(:,index(1));
norm = 1/sqrt(c'*S*c);
c = c*norm;

% anharmonic weighting
a=0;
% number of interations for scf
niter=10;
% stores the ground state energy for every iteration (for plotting)
evals = zeros(1, niter);

% Solving for ground state, iteratively
for iter=1:niter
    % anharmonic term
    avy4 = c'*G*c;
    % effective hamiltonian (mean field theory)
    heff = h + a*avy4*G;
    
    % Solve for new heff eigenvalues (energies) and eigenvectors
    % (bases coefficients)
    Sinvheff = Sinv*heff;
    [U D] = eig(Sinvheff);
    eigs = diag(D);
    [eigs index] = sort(eigs);
    
    % Get new ground state basis coeffs for next iteration
    c = U(:,index(1));
    norm = 1/sqrt(c'*S*c);
    c = c*norm;
    
    % save ground state energy for this iteration to observe convergence
    energy = 2*c'*h*c + a*(c'*G*c)^2;
    evals(iter) = energy;
end

plot(1:niter, evals)