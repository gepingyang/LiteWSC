%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo for the LiteWWSC algorithm, which is proposed in the %
% following paper:                                                  %
%                                                                   %
%LiteWSC:a Lightweight framework for Web-Scale Spectral Clustering %
%                                                                   %
% The code has been tested in Matlab R2018b                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ac_avg,nmi_avg, time_avg]=demoLiteWWSC()

dataname = {'./dataset/USPS.mat'}; %if you have many batch data, you should shuffer them and  put them in this cell in turn.                   
s =4000; %the number of data points that are used to generate prototypes and prototype graph.
p = 500; % the number of prototypes
r =5; %the number of nearest prototypes.
k =10; % the number of clusters.
seed.end = 10;% the number of times that the LiteWWSC runs.
seed.start = 1;
interval = seed.end - seed.start + 1;
ac_sum = 0;
nmi_sum = 0;
time_all = 0;
for i = seed.start : seed.end 
    rand('seed',i);    
     fprintf('Seed No: %d\n',i);
     tic; 
     [label_orig, label] = LiteWSC(dataname, k, s, p, r);
     time_once = toc;
     time_all = time_once + time_all;
     %%%%%%%%%%%%%%%%%%calculate clustering accuaracy%%%%%%%%%%%%%%%%%%%%%%
     label = bestMap(label_orig,label);
      nmi_result = nmi(label,label_orig);  
      ac_result = length(find(label_orig == label))/length(label);
      fprintf('nmi: %.2f%% +', nmi_result * 100);
      fprintf('ac: %.2f%% + ', ac_result * 100);
      fprintf("runtime: %.2f%% \n", time_once);
      nmi_sum = nmi_result + nmi_sum;
      ac_sum = ac_result + ac_sum;
end
nmi_avg = nmi_sum / interval;
ac_avg = ac_sum / interval;
time_avg = time_all / interval;
fprintf('**************************************************************\n');
fprintf('avg_nmi: %.2f%% + ', nmi_avg * 100);
fprintf('avg_ac: %.2f%% + ', ac_avg * 100);
fprintf('avg_runtime: %.2f s\n\n\n', time_avg);
fprintf('Algorithm Finished');
% end

     
