%% Setup VARIABLE parameters

Nside=7;
density = 0.7;

%% Setup DEPENDENT parameters
N = Nside^2;            % number of particles
L = sqrt(N/density);    % box dimension
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
drawconfig(r,N,L,Lmax)

% sets max magnitude of possible move per step
dx=1;

% initialize histogram
h.count=0;
h.range=[0 L/2];
h.binwidth=0.1;


%% MC Simulation

nsweeps=1000;
for sweep=1:nsweeps
    
  % For each sweep, attempt N random moves
  for step=1:N
    % Pick random particle and random move
    index = ceil(rand*N);
    rtrial = r(:,index) + (rand(2,1)-0.5)*dx;
    
    % Check to see if proposed move would result in overlap
    overlap=0;
    for i=1:N
      if (i ~= index && overlap==0)
        deltar = r(:,i)-rtrial;
        deltar = deltar - L*round(deltar/L);
        deltar2 = deltar'*deltar;
        if (deltar2 < 1)
          overlap=1;
        end
      end
    end
    
    % If no overlap, accept move
    if (overlap==0)
      r(:,index) = rtrial;
    end
  end
  
  % draw system every 100 sweeps
  if (mod(sweep,10)==0)
    drawconfig(r,N,L,Lmax)
  end
  
  % update histogram with all pairwise separation distances
  for i=1:N
    for j=i+1:N
      deltar = r(:,i)-r(:,j);
      deltar = deltar - L*round(deltar/L);
      deltar2 = deltar'*deltar;
      dr = sqrt(deltar2);
      h = histo(h,dr);
    end
  end
  
  
end
