clear
%clf

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

a=0.1;
niter=10;

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

energy
