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

chg = [2 1];

dist = 1.4632;
Rnuc = zeros(3,Nnuc);
Rnuc(1,2) = dist;
R = zeros(3,L,K);
for p=1:L
  R(1,p,2) = dist;
end


S = zeros(K,K);
T = zeros(K,K);
U1 = zeros(K,K);
U2 = zeros(K,K,K,K);

for mu=1:K
  for p=1:L
    al = alpha(p,mu);
    RA = R(:,p,mu)';
    for nu=1:K
      for q=1:L
        bet = alpha(q,nu);
        RB = R(:,q,nu)';
        
        S(mu,nu) = S(mu,nu) + d(p,mu)*d(q,nu)*overlap(al,bet,RA,RB);
        T(mu,nu) = T(mu,nu) + d(p,mu)*d(q,nu)*kinetic(al,bet,RA,RB);
        
        for n=1:Nnuc
          RC = Rnuc(:,n)';
          U1(mu,nu) = U1(mu,nu) + ...
            chg(n)*d(p,mu)*d(q,nu)*elec_nuc(al,bet,RA,RB,RC);
        end
        
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

h = T+U1;
Sinv = inv(S);
Sinvh = Sinv*h;
[U D] = eig(Sinvh);
eigs = diag(D);
[eigs index] = sort(eigs);

for i=1:N/2
  c(:,i) = U(:,index(i));
  norm = 1/sqrt(c(:,i)'*S*c(:,i));
  c(:,i) = c(:,i)*norm;
end
P = 2*c*c';

% clf
% xvals = -1:0.1:3;
% yvals = -2:0.1:2;
% sz = numel(xvals);
% psi = zeros(sz,sz);
% 
% for i=1:sz
%   x = xvals(i);
%   for j=1:sz
%     y = yvals(j);
%     for p=1:L
%       dx = x - R(1,p,1);
%       dy = y - R(2,p,1);
%       psi(j,i) = psi(j,i) + d(p,1)*exp(-alpha(p,1)*(dx^2+dy^2));
%     end
%   end
% end
% 
% contour(xvals,yvals,psi,10)
% axis equal
% hold on
% 
% plotcircle([0 0 0],0.2)
% plotcircle([dist 0 0],0.1)