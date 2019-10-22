x=[]; d=[]; chi=[]; allsigb=[]; allb=[];
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'Lgsv3.mat'))
    load(strcat('L',num2str(k),'.mat'))
    for kk=1:size(lgss,1)
        for kkk=1:size(lgss,2)
            if ~isempty(lgs{kk,kkk})
                lgs{kk,kkk}=lgs{kk,kkk}./data.fsample;
            pv=myBinomTest(length(find(lgs{kk,kkk}<0)),length(lgs{kk,kkk}),.5);
            allb = [allb lgs{kk,kkk}];
            if pv<=(.05/362)
                x=[x lgs{kk,kkk}./data.fsample];
                d=[d lgsm(kk,kkk)./data.fsample];
                chi=[chi; k kk kkk pv];
                allsigb = [allsigb lgs{kk,kkk}];
            end
            end
        end
    end
end
%%
figure
subplot(1,2,1)
hist(-allsigb,700), xlim([-.4 .4]), ylim([0 300]), box off
set(gca,'YTick',[0 150 300],'XTick',[-.7 -.35 0 .35 .7])
yL = get(gca,'YLim');
line([0 0],yL,'Color','r','LineWidth',2);
%line([median(allsigb) median(allsigb)],yL,'Color','b','LineWidth',2);
subplot(1,2,2)
hist(-d,10), xlim([-.7 .7]), box off
set(gca,'YTick',[0 5 10],'XTick',[-.7 -.35 0 .35 .7])
yL = get(gca,'YLim');
line([0 0],yL,'Color','r','LineWidth',2);
set(gcf,'PaperSize',[15 15])
% subplot(1,3,3)
% hist(allb,5000), xlim([-2.5 2.5]), box off
% yL = get(gca,'YLim');
% line([0 0],yL,'Color','r');