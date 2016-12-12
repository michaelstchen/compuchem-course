%% Set up parameters
len_OH1 = 1.809;                %OH #1 equilibrium bond len
len_OH2 = 1.809;                %OH #2 equilibrium bond len
ang_HOH = degtorad(104.52);     %HOH equilibrium bond angle
chg = [1 1 8];
N=10;                           %num electrons

%% Initial Calculations and Data Setup

Nnuc=3;     %num nuclei
K=7;        %num of basis functions
L=2;        %num of terms per basis

zeta = [1.24 7.66];     %zeta constants for H and O

alpha = zeros(L,K);     %basis fns alpha constats
d = zeros(L,K);         %basis fns term coeffs

alpha(:,1) = [0.151623 0.851819]*zeta(1)^2;     %1s H
alpha(:,2) = [0.151623 0.851819]*zeta(1)^2;     %1s H
alpha(:,3) = [0.151623 0.851819]*zeta(2)^2;     %1s O
alpha(:,4) = [0.493363 1.94523];                %2s O
alpha(:,5) = [0.9 0.9];                         %2px O
alpha(:,6) = [0.9 0.9];                         %2py O
alpha(:,7) = [0.9 0.9];                         %2pz O


d(:,1) = [0.164964 0.381381];   %1s H
d(:,2) = [0.164964 0.381381];   %1s H
d(:,3) = [0.164964 0.381381];   %1s O
d(:,4) = [0.168105 0.0241442];  %2s O
d(:,5) = [1 -1];                %2px O
d(:,6) = [1 -1];                %2py O
d(:,7) = [1 -1];                %2pz O
  
%Set nuclei coordinates
Rnuc = zeros(3,Nnuc);
Rnuc(:, 1) = [-len_OH1*sin(ang_HOH/2) -len_OH1*cos(ang_HOH/2) 0]; %Hydrogen1
Rnuc(:, 2) = [len_OH2*sin(ang_HOH/2) -len_OH2*cos(ang_HOH/2) 0];  %Hydrogen2
Rnuc(:, 3) = [0 0 0];                                             %Oxygen

% Set centers for basis function terms
dx = 0.1;
R = zeros(3,L,K);
% Hydrogen1 1s
R(:,1,1) = Rnuc(:, 1);
R(:,2,1) = Rnuc(:, 1);
% Hydrogen2 1s
R(:,1,2) = Rnuc(:, 2);
R(:,2,2) = Rnuc(:, 2);
% Oxygen 1s
R(:,1,3) = Rnuc(:, 3);
R(:,2,3) = Rnuc(:, 3);
% Oxygen 2s
R(:,1,4) = Rnuc(:, 3);
R(:,2,4) = Rnuc(:, 3);
% Oxygen 2px
R(:,1,5) = Rnuc(:, 3) + [dx 0 0]';
R(:,2,5) = Rnuc(:, 3) - [dx 0 0]';
% Oxygen 2py
R(:,1,6) = Rnuc(:, 3) + [0 dx 0]';
R(:,2,6) = Rnuc(:, 3) - [0 dx 0]';
% Oxygen 2pz
R(:,1,7) = Rnuc(:, 3) + [0 0 dx]';
R(:,2,7) = Rnuc(:, 3) - [0 0 dx]';


%% Initial Fock matrix guess
S = zeros(K,K);         %overlap matrix elements
T = zeros(K,K);         %kinetic matric elements
U1 = zeros(K,K);        %one-electron potential matrix elements
U2 = zeros(K,K,K,K);    %two-electron potential matrix elements

% Detrmining matrix elements. mu and nu are both bases.
% p and q are terms of the bases mu and nu, respectively
for mu=1:K
    for p=1:L
      al = alpha(p,mu);
      RA = R(:,p,mu)';
      for nu=1:K
        for q=1:L
          bet = alpha(q,nu);
          RB = R(:,q,nu)';

          %overlap matrix elements
          S(mu,nu) = S(mu,nu) + d(p,mu)*d(q,nu)*overlap(al,bet,RA,RB);
          %kinetic matrix elements
          T(mu,nu) = T(mu,nu) + d(p,mu)*d(q,nu)*kinetic(al,bet,RA,RB);

          %one-elec potential matrix elements (elec-nuc interactionss)
          for n=1:Nnuc
            RC = Rnuc(:,n)';
            U1(mu,nu) = U1(mu,nu) + ...
              chg(n)*d(p,mu)*d(q,nu)*elec_nuc(al,bet,RA,RB,RC);
          end

          %two-elec potential matrix elements (elec-elec interactions)
          % again we iterate over bases lambda and sigma with terms
          % s and t, respectively
          for lambda=1:K
            for s=1:L
              gam = alpha(s,lambda);
              RC = R(:,s,lambda)';
              for sigma=1:K
                for t=1:L
                  del = alpha(t,sigma);
                  RD = R(:,t,sigma)';

                  U2(mu,nu,lambda,sigma) = U2(mu,nu,lambda,sigma) + ...
                    d(p,mu)*d(q,nu)*d(s,lambda)*d(t,sigma)* ...
                    elec_elec(al,bet,gam,del,RA,RB,RC,RD);
                end
              end
            end
          end


        end
      end
    end
end

%Calculating energies and wavefn basis coeffs (one-electron)
h = T+U1;
Sinv = inv(S);
Sinvh = Sinv*h;
[U D] = eig(Sinvh);
eigs = diag(D);
[eigs index] = sort(eigs);

%% Iteratively improve Fock matrix
for i=1:N/2
    c(:,i) = U(:,index(i));
    norm = 1/sqrt(c(:,i)'*S*c(:,i));
    c(:,i) = c(:,i)*norm;
end
%our density matrix P
P = 2*c*c';

niter=5;
for iter=1:niter
    %calculating our Fock matrix F
    F = h;
    for mu=1:K
      for nu=1:K
        for lambda=1:K
          for sigma=1:K
            F(mu,nu) = F(mu,nu) + P(lambda,sigma)*...
              (U2(mu,nu,lambda,sigma) - 0.5*U2(mu,sigma,lambda,nu));
          end
        end
      end
    end

    %calculate fock eigenvalues and eigenvectorss
    SinvF = Sinv*F;
    [U D] = eig(SinvF);
    eigs = diag(D);
    [eigs index] = sort(eigs);

    %obtain our molecule's energies for this iteration
    energy = 0.5*trace(P*(h+F));

    %update density matrix for next iteration
    for i=1:N/2
      c(:,i) = U(:,index(i));
      norm = 1/sqrt(c(:,i)'*S*c(:,i));
      c(:,i) = c(:,i)*norm;
    end
    P = 2*c*c';

end

% add the nuclear repulsion terms
energy = energy + 8/len_OH1 + 8/len_OH2 ...
         + 1/(len_OH1*sin(ang_HOH/2)+len_OH2*sin(ang_HOH/2));
