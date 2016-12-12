% Used for Problem 2 parts (v) and (vi)
clear all;
close all;

M=30;               % Number of oscillators
E=30;               % Total amount of energy

n = zeros(1,M);     % Distribution of microstates

% Arbitarily initialize distrib
n(5)=E;

% nsteps = 100000;    % Number of simulation iterations

h.count=0;
h.range = [-0.5, E+0.5];
h.binwidth = 1;
h = histo(h,0);
p_actual = 0.25424;
eps = 0.01;
steps = 0;
while abs(h.hist(2)/h.count - p_actual) > eps
%for step=1:nsteps
  steps = steps + 1;
  i = ceil(rand*M);
  j = ceil(rand*M);
  
  if (n(i)>0)
    n(i) = n(i)-1;
    n(j) = n(j)+1;
  end
   
  for i = 1:length(n)
      h = histo(h,n(i));
  end
end

% bar(h.vals,h.hist/h.count)
% title('Counts of quanta in an oscillator');
% xlabel('n_i')
% ylabel('Count')
steps