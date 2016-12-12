function [ ground ] = CIOscfunc(N_orb, a)

% First use mean field to find configuration orbitals
alpha=2;
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

niter=100;
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

for i=1:size(U,2)
    norm = 1/sqrt(U(:,i)'*S*U(:,i));
    U(:,i) = U(:,i)*norm;
end


% Configuration Interaction estimate for 2D Non-linear Oscillator
N_conf = sum(1:N_orb-1);    %Number of CI configurations

% Enumerate the CI configurations
m = 1;
conf = zeros(2, N_conf);
for i = 1:N_orb
    for j = i+1:N_orb
        conf(1, m) = i;
        conf(2, m) = j;
        m = m + 1;
    end
end

% Determining the weights for the configurations
H_conf = zeros(N_conf, N_conf);
for li = 1:N_conf
    i = index(conf(1, li));
    j = index(conf(2, li));
    ci = U(:, i);
    cj = U(:, j);
    for ri = 1:N_conf
        k = index(conf(1, ri));
        l = index(conf(2, ri));
        ck = U(:, k);
        cl = U(:, l);
        
        H_conf(li, ri) = delta(i,k)*(cj'*h*cl) + delta(j,l)*(ci'*h*ck) ...
                - delta(i,l)*(cj'*h*ck) - delta(j,k)*(ci'*h*cl);
        H_conf(li, ri) = H_conf(li, ri) + a * ((ci'*G*ck)*(cj'*G*cl) ...
                                            - (cj'*G*ck)*(ci'*G*cl));
    end
end

[U_ci, D_ci] = eig(H_conf);
energs = sort(diag(D_ci));

ground = energs(1);

end

