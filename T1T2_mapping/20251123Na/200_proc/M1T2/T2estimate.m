clc;clear all;close all;
mask=BrikLoad('ROI+orig.BRIK');
data=[];BG=[];
raw='/home/xl872/Documents/raw_data/'; %raw data path to find the mathod file
dd=dir('*.an200epi+orig.BRIK'); %the anat data can not be name as xxxepixxx
time=[];
for d=1:length(dd)
    
    data0=squeeze(BrikLoad([dd(d).folder,'/',dd(d).name]));
     e=split(dd(d).name,'.');
    date=e{1};
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
    bg=data0([1:2,end-1:end],[1:2,end-1:end],:);
    BG=cat(2,BG,bg(:));
    data0=imgaussfilt3(data0,0.5);
    %soothes out the voxels using 3D gaussian filter with 0.5 sigma for the
    data=cat(4,data,data0);
    
    % concatenates the data with the dimenson of time(4)
end
[time,id]=sort(time);
data=data(:,:,:,id);
BG=BG(:,id);
%%
figure
imagesc(data(:,:,7,1)'); colormap("gray")
%%
figure
imagesc(mask(:,:,7)'); colormap("gray")
%%
figure;
subplot(2,1,1)
boxplot(BG);
xticklabels(num2str(time'))
ylabel('Image amplitude')
title('Backgroud value')
subplot(2,1,2)
scatter(time,std(BG),'filled');
ylabel('Backgroud STD')
xlabel('TE (s)')
xlim("tight")
fontsize(gcf,12,"points")
print(gcf,[pwd,'\Baseline.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\Baseline.eps'],'-depsc','-r300');
%% ROI mask
mask=BrikLoad('ROI+orig.BRIK');
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
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result.e
%
result1=result;tdata1=tdata;
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2fast=',num2str(result.b),'ms ','T2slow=',num2str(result.e),'ms'])
print(gcf,[pwd,'\aveT2_2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\aveT2_2estimate.eps'],'-depsc','-r300');
%

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,50, Inf],...
   'StartPoint',[tdata(1) 5 tdata(end)]);
ft = fittype('c+a*exp(-x/b)','options',fo); 
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result11=result;tdata11=tdata;
%
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2=',num2str(result.b),'ms'])
% legend('Location', 'southeast');
print(gcf,[pwd,'\aveT2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\aveT2estimate.eps'],'-depsc','-r300');

%% CSF mask
mask=BrikLoad('CSF+orig.BRIK');
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
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result.e
%
result2=result;tdata2=tdata;
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['CSFVT2fast=',num2str(result.b),'ms ','T2slow=',num2str(result.e),'ms'])
print(gcf,[pwd,'\CSFVaveT2_2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\CSFVaveT2_2estimate.eps'],'-depsc','-r300');
%

% tdata=adata/num;%not normalized
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,50, Inf],...
   'StartPoint',[tdata(1) 5 tdata(end)]);
ft = fittype('c+a*exp(-x/b)','options',fo); 
result = fit(time(1:end)', tdata(1:end),ft);
result.b
result22=result;tdata22=tdata;
%
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['CSFVT2=',num2str(result.b),'ms'])
% legend('Location', 'southeast');
print(gcf,[pwd,'\CSFVaveT2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\CSFVaveT2estimate.eps'],'-depsc','-r300');
%% compare ROI and CSF
figure(Position=[100 100 400 300]);
plot(result1,'r')
hold on
scatter(time(1:end)', tdata1(1:end),'r','.')

plot(result2,'b')
scatter(time(1:end)', tdata2(1:end),'b','.')

legend({'Cortex fitted','Cortex data','CSF fitted','CSF data'},'Box','off')
fontsize(gcf,10,"points")
ylabel('SNR')
xlabel('TE (ms)')
title({['CortexT2fast=',num2str(result1.b),'ms ','T2slow=',num2str(result1.e),'ms'] ...
    ['CSFT2fast=',num2str(result2.b),'ms ','T2slow=',num2str(result2.e),'ms']})
print(gcf,[pwd,'\aveT2_CSF&COTestimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\aveT2_CSF&COTestimate.eps'],'-depsc','-r300');
%

figure(Position=[100 100 400 300]);
plot(result11,'r')
hold on
scatter(time(1:end)', tdata1(1:end),'r','.')

plot(result22,'b')
scatter(time(1:end)', tdata2(1:end),'b','.')

legend({'Cortex fitted','Cortex data','CSF fitted','CSF data'},'Box','off')
fontsize(gcf,10,"points")
ylabel('SNR')
xlabel('TE (ms)')
title({['CortexT2=',num2str(result11.b),'ms '] ...
    ['CSFT2=',num2str(result22.b),'ms']})
print(gcf,[pwd,'\aveT2_CSF&COTestimate1.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\aveT2_CSF&COTestimate1.eps'],'-depsc','-r300');
%% each voxel
mask=BrikLoad('ROI+orig.BRIK');
maskCor=mask;
maskCSF=BrikLoad('ROI+orig.BRIK');
data=data(:,:,:,:);
T2=[];T2f=[];T2s=[];Tdata1=[];Tdata2=[];
for i=1:size(data,1)
    for j=1:size(data,2)
       for k=1:size(data,3)
            if mask(i,j,k)~=0
                % fo = fitoptions('Method','NonlinearLeastSquares',...
                %    'Lower',[0,0,0,0,0],...
                %    'Upper',[Inf,100, Inf,Inf,100],...
                %    'StartPoint',[tdata(1)/2 0.01 mdata/2 tdata(1)/100 10]);
                % ft = fittype('c+a*exp(-x/b)+d*exp(-x/e)','options',fo); 
                % % contains addition to the function to calculate T2
                % result = fit(time', squeeze(data(i,j,k,:)),ft);
                % T2f(i,j,k)=min(result.b,result.e);
                % T2s(i,j,k)=max(result.b,result.e);
                tdata=squeeze(data(i,j,k,:));%./data_std';
               if maskCor(i,j,k)==1
                    Tdata1=[Tdata1 tdata];
                elseif maskCSF(i,j,k)==1
                    Tdata2=[Tdata2 tdata];
                end
                fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0,0,5],...
                   'Upper',[Inf,5, Inf,Inf,100],...
                   'StartPoint',[tdata(1) 1 tdata(end) tdata(1) 5]);
                ft = fittype('c+a*exp(-x/b)+d*exp(-x/e)','options',fo); 
                result = fit(time(1:end)', tdata(1:end),ft);
                T2f(i,j,k)=min(result.b,result.e);
                T2s(i,j,k)=max(result.b,result.e);

                
                fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0],...
                   'Upper',[Inf,100, Inf],...
                   'StartPoint',[tdata(1) 5 tdata(end)]);
                ft = fittype('c+a*exp(-x/b)','options',fo); 
                result = fit(time(1:end)', tdata(1:end),ft);
                T2(i,j,k)=result.b;
            else
                T2f(i,j,k)=0;
                T2s(i,j,k)=0;
                T2(i,j,k)=0;
            end

       end
    end
end

%%
anat=squeeze(BrikLoad('ave_temp+orig.BRIK'));

% maskt=BrikLoad('COPY_20250101Na.13.an200epi.resamp+orig');
% maskt(:,80:end,:)=0;
mask=BrikLoad('ROI+orig.BRIK');
maskt=mask;
for sl=5:7%:size(data,3)
figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
B=anat(:,:,sl)';
B=repmat(B,1,1,3)/max(B(:));
hB = image(B);%axis xy
axis image off;
hold on;
F=T2f(:,:,sl)';
% F=imresize(F,[size(B,1),size(B,2)]);
hF = imagesc(F);%axis xy
colormap('jet')
clim([0 2]);
%
alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2*fast (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2festimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2festimate',num2str(sl),'.eps'],'-depsc','-r300');


figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
% B=anat(:,:,sl)';
% B=repmat(B,1,1,3)/max(B(:))*1.5;
hB = image(B);%axis xy
axis image off;
hold on;
F=T2s(:,:,sl)';
% F=imresize(F,[size(B,1),size(B,2)],"bilinear");
hF = imagesc(F);%axis xy
colormap('jet')
clim([0 15]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2*slow (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2sestimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2sestimate',num2str(sl),'.eps'],'-depsc','-r300');

figure(Position=[100 100 400 300]);
% B=data(:,:,sl,1)';
% B=anat(:,:,sl)';
% B=repmat(B,1,1,3)/max(B(:))*1.5;
hB = image(B);%axis xy
axis image off;
hold on;
F=T2(:,:,sl)';
% F=imresize(F,[size(B,1),size(B,2)],"bilinear");
hF = imagesc(F);%axis xy
colormap('jet')
clim([0 15]);
% alphadata =0.25*maskt(:,:,sl)';% 
set(hF,'AlphaData',alphadata);
h = colorbar;
ylabel(h, 'T2* (ms)');
fontsize(gcf,12,"points")

print(gcf,[pwd,'/T2estimate',num2str(sl),'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'/T2estimate',num2str(sl),'.eps'],'-depsc','-r300');
end

%%
figure;


 
 fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0,0,0,5],...
                   'Upper',[Inf,5, Inf,Inf,100],...
                   'StartPoint',[tdata(1) 0.1 tdata(end) tdata(1) 5]);
ft = fittype('c+a*exp(-x/b)+d*exp(-x/e)','options',fo); 
result1 = fit(repmat(time(1:end)',[size(Tdata1,2),1]), Tdata1(:),ft);

result2 = fit(repmat(time(1:end)',[size(Tdata2,2),1]), Tdata2(:),ft);

plot(result1)
hold on
plot(result2,'b')
scatter(time,Tdata1,'r','filled');

scatter(time,Tdata2,'b','filled');
title({['Cortex T2_f=',num2str(result1.b),'ms T2_s=',num2str(result1.e),'ms'] ...
    ['CSF T2_f=',num2str(result2.b),'ms T2_s=',num2str(result2.e),'ms'] })
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
legend({'Cortex fitted result' 'CSF fitted result'});
print(gcf,[pwd,'\allT2estimate.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\allT2estimate.eps'],'-depsc','-r300');


%%
sl=9;

figure;imagesc(T2(:,:,sl)')
tdata=squeeze(data(9,10,sl,:));
fo = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0,0,0],...
   'Upper',[Inf,50, Inf],...
   'StartPoint',[tdata(1) 5 tdata(end)]);
ft = fittype('c+a*exp(-x/b)','options',fo); 
result = fit(time(1:end)', tdata(1:end),ft);
result.b
%
figure(Position=[100 100 400 300]);
plot(result,time(1:end)', tdata(1:end))
fontsize(gcf,12,"points")
ylabel('SNR')
xlabel('TE (ms)')
title(['T2=',num2str(result.b),'ms'])