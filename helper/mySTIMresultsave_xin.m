function mySTIMresultsave(k,bk1,bk2,stimtime,TR,xlim1,xlim2,fs,aveTrial,MC,Nt,MD,I,timepoint,timepoint2,name)
s=k;
ss=1.2*s;%step length
step=round(1*s*3);%half window
Ca=[];Na=[];ii=1;Ca_ave=[];
% colors = orderedcolors("gem");
for i=step:ss:(length(I)-step)
    %based on GLU
  Ca(ii)=mean(MC(I(i-step+1:i+step))); 
  Na(ii,:)=mean(Nt(I(i-step+1:i+step),:),1);
  Ca_ave(ii,:)=mean(MD(:,I(i-step+1:i+step)),2);
  % based on order
  % Ca(ii)=mean(MC(i-step+1:i+step));  
  % Na(ii,:)=mean(Nt((i-step+1:i+step),:),1);
  ii=ii+1;
end
Naout=100*(mean(Na(:,timepoint:timepoint2),2)-mean(Na(:,bk1:bk2),2))./mean(Na(:,bk1:bk2),2);
Caout=zscore(Ca)';
figure(Position=[10 50 1000 300]);
subplot(1,3,1)

scatter(Caout,Naout,'k','filled')

mdl = fitlm(Caout,Naout);
 hold on
plot(Caout,mdl.Fitted,LineWidth=1,Color='r')
xlabel('Z-scaled Glu change')
ylabel(['NARS-fMRI change(%)'])
 % title({['Mean of ',num2str(2*step*aveTrial),' trials'],['k=',num2str(mdl.Coefficients.Estimate(2))],['p=',num2str(mdl.Coefficients.pValue(2))]})
xlim('tight')
subplot(1,3,2)
n=ii-1;


colors = jet(n+1); % 你可以试试 parula, hot, cool, etc.
t=1:size(MD,1);
t=t/fs;
for i = 1:n  
     h = plot(t,Ca_ave(i,:)-Ca_ave(i,1), 'Color', [colors(i+1,:) 1],'LineWidth',1.5);  % 添加透明度通道
     hold on
end
xlabel('Time (s)')
ylabel('Normalized Glu (%)')
set(gca,'box','off')
subplot(1,3,3)
t=1:size(Na,2);
t=t-stimtime;
t=t*TR;
for i = 1:n  
     h = plot(t,Na(i,:), 'Color', [colors(i+1,:) 1],'LineWidth',1.5);  % 添加透明度通道
     hold on
end
xlabel('Time (s)')
ylabel('NARS-fMRI (%)')
set(gca,'box','off')
fontsize(gcf,18,"points")
xlim(TR*[xlim1-stimtime,xlim2-stimtime])
print(gcf,[pwd,'\STIM',name,'.jpg'],'-djpeg','-r300');
print(gcf,[pwd,'\STIM',name,'.eps'],'-depsc','-r300');
save(['STIM',name,'.mat'],"Caout","Naout")
T=table(Caout,Naout);
writetable(T,['STIM',name,'.xlsx'],'Sheet',1);
end