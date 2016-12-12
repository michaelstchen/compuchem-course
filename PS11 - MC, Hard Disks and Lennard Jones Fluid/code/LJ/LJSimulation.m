%% Setup parameters
density = 0.75;
T = 0.45;
beta = 1 / T;

Nside=5;
N = Nside^2;            % number of particles
L = max(sqrt(N/density), 7);    % box dimension
r = zeros(2,N);         % particle positions

% Initialize particle configurations
spacing = L/Nside;
index=0;
for i=1:Nside
  for j=1:Nside
    index = index+1;
    x = i*spacing + spacing/2;
    y = j*spacing + spacing/2;
    r(:,index) = [x y]';
  end
end

Lmax=10;
% drawconfig(r,N,L,Lmax);

% initial energy
U = LJPotentialTotal(N, L, r);

% sets max magnitude of possible move per step
dx=1;

% initialize histogram
hr.count=0;
hr.range=[0, L/2];
hr.binwidth=0.1;

hu.count = 0;
hu.range = [-50, -30];
hu.binwidth=0.1;

%% MC Simulation

nsweeps=10000;
for sweep=1:nsweeps
    
  % For each sweep, attempt N random moves
  for step=1:N
    % Pick random particle and random move
    index = ceil(rand*N);
    rorig = r(:, index);
    rtrial = rorig + (rand(2,1)-0.5)*dx;
    
    % Calculate change in energy
    dU = 0;
    overlap=0;
    for i=1:N
      if (i ~= index && overlap==0)
        deltar_orig = r(:,i)-rorig;
        deltar_orig = deltar_orig - L*round(deltar_orig/L);
        dr_orig = sqrt(deltar_orig'*deltar_orig);
          
        deltar_trial = r(:,i)-rtrial;
        deltar_trial = deltar_trial - L*round(deltar_trial/L);
        dr_trial = sqrt(deltar_trial'*deltar_trial);
        if (dr_trial < 1)
          overlap=1;
        end
        
        dU = dU + LJPotential(dr_trial, 2.5) - LJPotential(dr_orig, 2.5);
      end
    end
    
    % If no overlap, accept move
    if (overlap==0 && rand < exp(-beta*dU))
      r(:,index) = rtrial;
      U = U + dU;
    end
  end
  
  % draw system every 100 sweeps
  if (mod(sweep,10)==0)
    drawconfig(r,N,L,Lmax)
  end
  
  % update histogram with all pairwise separation distances
  for i=1:N
    for j=i+1:N
      deltar_trial = r(:,i)-r(:,j);
      deltar_trial = deltar_trial - L*round(deltar_trial/L);
      deltar2 = deltar_trial'*deltar_trial;
      dr = sqrt(deltar2);
      hr = histo(hr, dr);
    end
  end
  
  hu = histo(hu, U);
  
end
