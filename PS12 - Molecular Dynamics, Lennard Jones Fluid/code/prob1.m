%% Problem 1 part(ii)
clear all;
close all;
clc;

N = 25;
L = 7;

r = zeros(2,N);

spacing = L/sqrt(N);
index=0;
for i=1:sqrt(N)
  for j=1:sqrt(N)
    x = i*spacing + spacing/2;
    y = j*spacing + spacing/2;
    
    index = index+1;
    r(:,index) = [x y]';
  end
end

drawconfig(r,N,L,L+0.5);

[F, U] = LJForcePot(r, 2.5, L);
fprintf('Force for one particle: [%.3f, %.3f]\n', F(1, 1), F(2, 1));

%% Problem 1 part (iii)
clear all;
close all;
clc;

N = 25;
L = 7;

r = zeros(2,N);

spacing = L/sqrt(N);
index=0;
for i=1:sqrt(N)
  for j=1:sqrt(N)
    x = i*spacing + spacing/2;
    y = j*spacing + spacing/2;
    
    index = index+1;
    r(:,index) = [x y]';
  end
end

r(1, 1) = r(1, 1) + 0.1;

drawconfig(r,N,L,L+0.5);

[F, U] = LJForcePot(r, 2.5, L);
fprintf('Force for shifted: [%.3f, %.3f]\n', F(1, 1), F(2, 1));
fprintf('Force for neighbor: [%.3f, %.3f]\n', F(1, 2), F(2, 2));
