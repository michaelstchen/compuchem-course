% %% Problem 1 part (vi)
% clear all;
% close all;
% clc;
% 
% plist = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
% avden_list = zeros(1, length(plist));
% 
% for i = 1:length(plist)
%     clearvars -except plist avden_list i
%     
%     pressure = plist(i);
%     density = 0.7;
%     Nside = 5;
%     
%     echodemo('hard_disks.m', 2);
%     echodemo('hard_disks.m', 3);
%     
%     avden_list(i) = avden;
% end
% 
% plot(plist, avden_list)
% xlabel('pressure')
% ylabel('<\rho>')
close
d = [0.8727,0.8901,0.895,0.8979,0.9068,0.9114,0.9148,0.92,0.925,0.927];
p = [9,10,10.5,11,11.25,11.5,12,13,14,15];
plot(d, p)