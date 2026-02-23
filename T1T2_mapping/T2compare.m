%%
clc;clear all;close all
load T2all_cortexfix2.mat;
T2fix=T2all;
load T2all_cortex2.mat;
data3=[T2fix.cortex_r2; T2all.cortex_r2];
figure(Position=[100 100 400 250]);
boxPlot(data3', ...
    'boxLabels', {'Fixed' 'Not Fixed'}, ...
    'plotPoints', true,'pointSize', 100);
ylabel('R2')
% ylim([0.98 1])
fontsize(gcf,20,"points")
print(gcf,['CortexT2_compare2estimate_r2.jpg'],'-djpeg','-r300');
