clc;clear all;close all;

raw='/home/xl872/Documents/raw_data/';
dall=dir('/home/xy079/T1T2/upload/2025*/200_proc/*T1');

T1all=struct;
figure(Position=[100 100 400 300]);
for da=1:length(dall)
    % close all
    dd=dir([dall(da).folder,'/',dall(da).name,'/*.an200epi+orig.BRIK']); %the anat data can not be name as xxxepixxx
    time=[];data=[];BG=[];
    for d=1:length(dd)
        %raw data path to find the mathod file
        data0=squeeze(BrikLoad([dd(d).folder,'/',dd(d).name]));
         e=split(dd(d).name,'.');
        date=e{1};
        T1all(da).date=[date,dall(da).name];
        enm=e{2};
        filename = ([raw,date,'/',enm,'/method']);
        fid = fopen(filename,'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        lines=lines{:};
        lineidx=find(strncmp(lines,'##$PVM_RepetitionTime=',16));
        B0=split([lines{lineidx}],'=');
        B=str2num(B0{2});
        time(d)=B;
        
        % loads a series of data with a number to string conversion
        
        %filter size of 3 (2*ceil(2*sigma)+1
        bg=data0(:,[1:1,end-1:end],:);
        
        BG=cat(2,BG,bg(:));
        data0=imgaussfilt3(data0,0.5);
        %soothes out the voxels using 3D gaussian filter with 0.5 sigma for the
        data=cat(4,data,data0);
        
        % concatenates the data with the dimenson of time(4)
    end
    [time,id]=sort(time);
    data=data(:,:,:,id);
    BG=BG(:,id);

    %% Cortex
    mask=BrikLoad([dd(d).folder,'/ROI+orig.BRIK']);
    maskCSF=BrikLoad([dd(d).folder,'/CSF+orig.BRIK']);
    mask(maskCSF==1)=0;
    data_std=std(BG);
    endn=length(time);
    %amount of time points we want to use
    
    data1=data(:,:,:,1:endn);
    data1=data1.*repmat(mask,1,1,1,size(data1,4));
    tdata=squeeze(sum(data1,[1,2,3]))/sum(mask,"all")./std(BG)';
    
    % tdata=adata/num;%not normalized
    fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,100, Inf],...
   'StartPoint',[tdata(1) 50 tdata(end)]);
    ft = fittype('c+a*(1-exp(-x/b))','options',fo); 

    [result,gof] = fit(time(1:end)', tdata(1:end),ft);
    result.b
    result1=result;tdata1=tdata;
    %
    T1all(da).cortex_r2=gof.rsquare;
    plot(result,time(1:end)', tdata(1:end))
    hold on
    fontsize(gcf,12,"points")
    ylabel('SNR')
    xlabel('TR (ms)')
    % title(['T1=',num2str(result.b),'ms'])
    legend('off');
    T1all(da).cortex=result.b;
    
  
end
save T1all T1all
xlim([0 275])
fontsize(gcf,20,"points")
print(gcf,['CortexT1estimateAll.jpg'],'-djpeg','-r300');
print(gcf,['CortexT1estimateAll.eps'],'-depsc','-r300');
T = struct2table(T1all);
writetable(T, 'T1all.csv');
%%
data1=[T1all.cortex]';
data2=[T1all.cortex_r2]';

figure(Position=[100 100 500 300]);
subplot(1,2,1)
boxPlot(data1, ...
    'boxLabels', {'Cortex'}, ...
    'plotPoints', true,'pointSize', 100);
ylabel('T1 (ms)')

subplot(1,2,2)
boxPlot(data2, ...
    'boxLabels', {'Cortex'}, ...
    'plotPoints', true,'pointSize', 100);
ylabel('R2')
fontsize(gcf,20,"points")
print(gcf,['CortexT1.jpg'],'-djpeg','-r300');
