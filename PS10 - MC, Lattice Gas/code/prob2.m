% %% Problem 2 part(ii)
% clear all;
% close all;
% clc;
% 
% L = 10;
% M = L*L;
% 
% mu = -2;
% T = 0.8;
% beta = 1/T;
% 
% nsteps = 10^6;
% 
% h.count=0;
% h.range=[-0.5 M+0.5];
% h.binwidth=1; 
% 
% echodemo('latticegas_interacting.m', 2);
% 
% plot(h.vals, h.hist/h.count)
% title(sprintf('T = %.1f', T));
% xlabel('N');
% ylabel('P(N)');


% % %% Problem 2 part(iii) + (iv)
% % clear all;
% % close all;
% % clc;
% % 
% % L = 10;
% % M = L*L;
% % 
% % mu = -2;
% % T_list = 0.4:0.1:1.0;
% % nsteps = 10^6;
% % 
% % for i = 1:length(T_list)
% %     T = T_list(i);
% %     beta = 1 / T;
% %     
% %     h.count=0;
% %     h.range=[-0.5 M+0.5];
% %     h.binwidth=1; 
% %     h.numbins=101;
% %     h.hist = [];
% %     h.vals = [];   
% %     
% %     echodemo('latticegas_interacting.m', 2);
% % 
% %     subplot(3, 3, i);
% %     plot(h.vals, h.hist/h.count)
% %     title(sprintf('T = %.1f', T));
% %     xlabel('N');
% %     xlim([0 25])
% %     ylabel('P(N)');
% %     
% %     probs = h.hist / h.count;
% %     norm_probs = probs / max(probs);
% %     plot(h.vals, log(norm_probs) / M);
% %     title(sprintf('T = %.1f', T));
% %     xlabel('N');
% %     xlim([0 25])
% %     ylabel('ln(P(N)/P(N*)) / M');
% % end

%% Problem 2 part(vii)

L = 10;
M = L*L;

mu = -2;
T = 0.4:0.1:1.0;
beta=1./T;
nreps = length(T);

n = zeros(L,L,nreps);
N=zeros(1,nreps);
E=zeros(1,nreps);
positions = zeros(2,M,nreps);

nsteps=10^6;
    
hist_list = struct();
for i = 1:length(T)
    h_list(i).count = 0;
    h_list(i).range=[-0.5 M+0.5];
    h_list(i).binwidth=1;
    h_list(i).numbins=101;
    h_list(i).hist = [];
    h_list(i).vals = [];    
end

echodemo('latticegas_parallel.m', 2);

for i = 1:length(T)
    subplot(3, 3, i)
    h = h_list(i);
    plot(h.vals, h.hist/h.count);
    title(sprintf('T = %.1f', T(i)));
    xlabel('N');
    xlim([0 M])
    ylabel('P(N)');
end
