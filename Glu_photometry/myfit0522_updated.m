clc;clear all;close all;
%% loas Glu
bionameall={'23Na', 'epi'};
biochannel=7;triggerchannel=4;ntrigger=225;
cutlength=0.2;%0.2s
[AD,trignum,fs]=myloadBio(bionameall,biochannel,triggerchannel,ntrigger,cutlength);

%% load mask&data 
mask=BrikLoad([pwd,'/23Na_ave_Glu/ROI_LFP_merge+orig.BRIK']);

date='05222025';
enum=[27 38; 41 69; 72 91; 95 134; 137 185; 189 192; 194 225; 228 246;];
Nt=myloadNa(date,mask,enum);

%%
save Na0522_channel4_new2.mat Nt AD fs

%% load Glu
clc;clear all;close all;
ntrigger=225;
cutlength=0.2;%0.2s
name='Na0522_channel4_new2';
load([name,'.mat'])
timepoint=11;
timepoint2=12;
GlustimNum=225;
aveTrial=1;
%% Mean Glu Mean Na

AD1=cell2mat(AD')';
AD1=reshape(AD1,size(AD1,1),ntrigger,size(AD1,2)/ntrigger);
MD=squeeze(median(AD1(:,1:GlustimNum,:),2));
% data average
CUT=floor(size(MD,2)/aveTrial)*aveTrial;
MD=MD(:,1:CUT);
Nt=Nt(1:CUT,:)';
MD=squeeze(mean(reshape(MD,[size(MD,1),CUT/aveTrial,aveTrial]),3));
Nt=squeeze(mean(reshape(Nt,[size(Nt,1),CUT/aveTrial,aveTrial]),3))';
%% ORDER (mean)
% MC=max(MD,[],1)-MD(1,:);
MC=mean(MD(0.03*fs:0.045*fs,:),1)-MD(1,:);
step=1;%half window
ss=1;%step length
Ca=[];Na=[];ii=1;
for i=step:ss:(length(MC)-step)
  Ca(ii)=mean(MC(i-step+1:i+step));  
  Na(ii,:)=mean(Nt((i-step+1:i+step),:),1);
  ii=ii+1;
end

tNa=100*(mean(Na(:,timepoint:timepoint2),2)-mean(Na(:,5:29),2))./mean(Na(:,5:29),2);


%% different step (mean)

% MC=max(MD,[],1)-MD(1,:);
MC=mean(MD(0.03*fs:0.045*fs,:),1)-MD(1,:);

[B,I]=sort(MC);



figure(Position=[10 50 1500 600]);
k=4;
for s=1:k
    ss=8*s;%step length
    step=1*s*12;%half window
    Ca=[];Na=[];ii=1;
    for i=step:ss:(length(B)-step)
        %based on GLU
      Ca(ii)=mean(MC(I(i-step+1:i+step))); 
      Na(ii,:)=mean(Nt(I(i-step+1:i+step),:),1);

      % based on order
      % Ca(ii)=mean(MC(i-step+1:i+step));  
      % Na(ii,:)=mean(Nt((i-step+1:i+step),:),1);
      ii=ii+1;
    end

    subplot(2,k,s)
    scatter(Ca,mean(Na(:,timepoint:timepoint2),2),'filled')
    xlabel('Glu (%)')
    ylabel(['Na (',num2str(timepoint),'-',num2str(timepoint2),')'])
    xlim('tight')
    
    fontsize(gcf,16,"points")
    mdl = fitlm(Ca,mean(Na(:,timepoint:timepoint2),2));
    % Gall(i1,i2,i3,1)=mdl.Coefficients.Estimate(2);
    % mdl.Coefficients.pValue(2)
    % GR2all(i1,i2,i3,1)=mdl.Rsquared.Ordinary;
    hold on
    plot(Ca,mdl.Fitted,LineWidth=1)
    title({['Mean of ',num2str(2*step*aveTrial),' trials'],['k=',num2str(mdl.Coefficients.Estimate(2))],['p=',num2str(mdl.Coefficients.pValue(2))]})
    
    subplot(2,k,k+s)
    scatter(Ca,100*(mean(Na(:,timepoint:timepoint2),2)-mean(Na(:,5:29),2))./mean(Na(:,5:29),2),'filled')
    mdl = fitlm(Ca,100*(mean(Na(:,timepoint:timepoint2),2)-mean(Na(:,5:29),2))./mean(Na(:,5:29),2));
     hold on
    plot(Ca,mdl.Fitted,LineWidth=1)
    xlabel('Glu (%)')
    ylabel(['Na (',num2str(timepoint),'-',num2str(timepoint2),') (f-f0)/f0 (%)'])
    title({['k=',num2str(mdl.Coefficients.Estimate(2))],['p=',num2str(mdl.Coefficients.pValue(2))]})
   
    xlim('tight')
    
end
fontsize(gcf,16,"points")
print(gcf,[pwd,'\Glu-Na',name,'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\Glu-Na',name,'.eps'],'-depsc','-r300');
%%
addpath('/home/xy079/Xin_proc/')
k=20;
bk1=5;bk2=29;TR=0.01;xlim1=5;xlim2=25;sttime=10;
 mySTIMresultsave_xin(k,bk1,bk2,sttime,TR,xlim1,xlim2,fs,aveTrial,MC,Nt,MD,I,timepoint,timepoint2,name)
%%
n = size(MD,2); % 总共要画多少条曲线

% 生成 colormap（例如从蓝到红）
colors = jet(n); % 你可以试试 parula, hot, cool, etc.
t=1:size(MD,1);
t=t/fs;
figure(Position=[100 100 550 300]);
hold on;
for i = 1:n  
     h = plot(t,MD(:,i)-MD(1,i), 'Color', [colors(i,:) 0.3]);  % 添加透明度通道
end
hold off;
colormap(jet(n));
cb=colorbar; % 可选，显示颜色条
caxis([1 n])  

cb.Ticks = 1:20:n;               % 设置刻度
cb.TickLabels = string(1:20:n); % 设置刻度标签（可选）
xlabel('Time (s)')
ylabel('Normalized Glu (%)')
cb.Label.String = 'Trial index';
xlim("tight")
fontsize(gcf,16,"points")
print(gcf,[pwd,'\AveragedGlu_all',name,'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\AveragedGlu_all',name,'.eps'],'-depsc','-r300');
%%
function Nt=myloadNa(date,mask,enum)
    Nt=[];tn=1;
    for i=1:size(enum,1)
        for d=enum(i,1):enum(i,2)
            dd=dir([pwd,'/23Na_ave_Glu/',date,'.',num2str(d),'.norm+orig.BRIK']);
            if length(dd)==1
            V=BrikLoad([dd.folder,'/',dd.name]);
            
            Nt(tn,:)=squeeze(sum(V.*repmat(mask,1,1,1,size(V,4)),[1,2,3]))/sum(mask,[1 2 3]);
            tn=tn+1;
            end
        end
    end
end
function [AD,trignum,fs]=myloadBio(bionameall,biochannel,triggerchannel,ntrigger,cutlength)
    AD={};trignum=[];tn=1;
    %just for 09192025
    delay=0;
    % if triggerchannel==4
    %     triggerchannel=2;delay=-0.4;
    % end
    for bio=bionameall
        bioname=bio{:};
        dd=dir([bioname,'*.mat']);
        ADt={};
        for d=1:length(dd)
            load(dd(d).name)
            fs=1000/isi;
            ad=[];
            cal=-data(:,biochannel);
            
            cal=(cal-mean(cal(1:fs)))/mean(cal(1:fs))*100;
            cal=detrend(cal,3);
            trigger=find(data(:,triggerchannel)>1);
            dt=diff([0; trigger]);
            trigger=trigger(find(dt~=1));
            trigger=trigger+delay*fs;%just for 09192025
            trignum(tn)=length(trigger)/ntrigger;tn=tn+1;
            for i=1:length(trigger)
                ad(i,:)=cal(trigger(i):trigger(i)+fs*cutlength-1);
            end
            ADt{d}=ad;
        end
        AD={AD{:},ADt{:}};
    end
end