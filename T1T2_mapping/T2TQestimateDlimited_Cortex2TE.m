clc;clear all;close all;
raw='/home/xl872/Documents/raw_data/';
dall=dir('/home/xy079/T1T2/upload/2025*/200_proc/*T2*');
output='/home/xy079/T1T2/upload/T2/';


T2all=struct;
figure(Position=[100 100 400 300]);
for da=1:length(dall)
    % close all
    dd=dir([dall(da).folder,'/',dall(da).name,'/2025*.an200epi+orig.BRIK']); %the anat data can not be name as xxxepixxx
    time=[];data=[];BG=[];
    for d=1:length(dd)
        %raw data path to find the mathod file
        data0=squeeze(BrikLoad([dd(d).folder,'/',dd(d).name]));
         e=split(dd(d).name,'.');
        date=e{1};
        T2all(da).date=[date,dall(da).name];
        enm=e{2};
        filename = ([raw,date,'/',enm,'/method']);
        fid = fopen(filename,'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        lines=lines{:};
        lineidx=find(strncmp(lines,'##$PVM_EchoTime=',15));
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
   'Lower',[0,0,0,0,5],...
   'Upper',[Inf,5, Inf,1,100],...
   'StartPoint',[tdata(1) 1 tdata(end) tdata(1) 5]);
    ft = fittype('c+a*((1-d)*exp(-x/b)+d*exp(-x/e))','options',fo); 
    [result,gof] = fit(time(1:end)', tdata(1:end),ft);
    result.b
    result.e
    T2all(da).cortex_r2=gof.rsquare;
    result1=result;tdata1=tdata;
    %
    % figure(Position=[100 100 400 300]);
    plot(result,time(1:end)', tdata(1:end))
    hold on
    fontsize(gcf,12,"points")
    ylabel('SNR')
    xlabel('TE (ms)')
    T2all(da).cortex_fast=result.b;
    T2all(da).cortex_slow=result.e;
    T2all(da).cortex_Weightslow=(result.d);
    T2all(da).cortex_Weightfast=(1-result.d);

  
      
end
save T2all_cortex2 T2all
legend('off')
fontsize(gcf,20,"points")
print(gcf,['CortexT2_2estimate.jpg'],'-djpeg','-r300');
T = struct2table(T2all);
writetable(T, 'T2all.csv');

%%
data1=[T2all.cortex_fast;T2all.cortex_slow]';
data2=[T2all.cortex_Weightfast;T2all.cortex_Weightslow]';

figure(Position=[100 100 600 300]);
subplot(1,2,1)
boxPlot(data1, ...
    'boxLabels', {'fast', 'slow'}, ...
    'plotPoints', true,'pointSize', 100);
ylabel('T2 (ms)')

subplot(1,2,2)
boxPlot(data2, ...
    'boxLabels', {'fast', 'slow'}, ...
    'plotPoints', true,'pointSize', 100);
ylabel('T2 weight')
fontsize(gcf,20,"points")
print(gcf,['CortexT2_2estimate_T2.jpg'],'-djpeg','-r300');

data3=[T2all.cortex_r2];
figure(Position=[100 100 200 300]);
boxPlot(data3', ...
    'boxLabels', {'Not Fixed'}, ...
    'plotPoints', true,'pointSize', 100);
ylabel('R2')
ylim([0.98 1])
fontsize(gcf,20,"points")
print(gcf,['CortexT2_2estimate_r2.jpg'],'-djpeg','-r300');
