%% Setup VARIABLE Parameters
Nside=5;
pressure = 11.25;
density = 0.90;

%% Setup DEPENDENT Parameters
N = Nside^2;            % number of particles

L = sqrt(N/density);    % box dimension
A = L^2;                % box area

r = zeros(2,N);         % particle positions

% Initialize particle configurations
spacing = L/Nside;
index=0;
for i=1:Nside
  for j=1:Nside
    index = index+1;
    x = i*spacing;
    y = j*spacing;
    r(:,index) = [x y]';
  end
end

Lmax=10;
% drawconfig(r,N,L,Lmax)

dx=1;               % determines max magnitude of possible move per step
b = 0.5;            % determines max magnitude of area change
nacc = 0;           % number of accepted moves
avden=0;            % keeps track of system's average density
nsweeps=30000;      % number of sweeps the simulation will run for

%% MC Simulation
for sweep=1:nsweeps
  
  % propose random change in box area / scale particle positions
  Atrial = A + (rand-0.5)*b;
  Ltrial = sqrt(Atrial);
  gamma = sqrt(Atrial/A);
  rtrial = gamma*r;
  
  % check new scaled positions do not overlap
  overlap=0;
  for i=1:N
    for j=i+1:N
      if (overlap==0)
        deltar = rtrial(:,i)-rtrial(:,j);
        deltar = deltar - Ltrial*round(deltar/Ltrial);
        deltar2 = deltar'*deltar;
        if (deltar2 < 1)
          overlap=1;
        end
      end
    end
  end
  
  % If no overlap and distribution obeyed, accept move
  if (overlap==0 && rand < ...
      exp(-pressure*(Atrial-A) + 2*N*log(gamma) ))
    r = rtrial;
    A = Atrial;
    L = Ltrial;
    nacc = nacc+1;
  end
    
  % For each sweep, attempt N random moves
  for step=1:N
    % Pick random particle and random move
    index = ceil(rand*N);
    rtrial = r(:,index) + (rand(2,1)-0.5)*dx;
    
    % check to see if proposed move would result in overlap
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
  if (mod(sweep,100)==0)
    drawconfig(r,N,L,Lmax)
  end
  
  % update average density
  if (sweep > 10000)
    avden = avden + N/A;
  end
end

avden = avden/(nsweeps-10000)      % calculate average density
facc = nacc/nsweeps;        % calculate acceptance frequency
