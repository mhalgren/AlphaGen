clear
s_ind=1;
%1.6.2 .2.61
cs=[0,0,0.562500000000000;0,0.875000000000000,1;1,0.812500000000000,0;0.500000000000000,0,0;];
figure
for k=[1 2 3 10 13 14 15 51 52]

    load(strcat('L',num2str(k),'Resultsv5.mat'),'freq','freq_ax','lbl')
    load(strcat('L',num2str(k),'ThChs.mat'))
    [~,fll]=min(abs(freq.freq-2));
    [~,ful]=min(abs(freq.freq-45));
    pow=freq.fourierspctrm(:,1:length(lbl),:); pow=pow.*conj(pow); 
    tmp=pow(:,:,fll:ful); %tmp=tmp(:);
    for kk=1:size(pow,2)
        tmp2=pow(:,kk,fll:ful);
        pow(:,kk,:)=pow(:,kk,:)./repmat(mean(tmp2(:)),[size(pow,1) 1 size(pow,3)]);
    end
    c_ind=1;
    for kk=ThChs
        subplot(3,3,s_ind)
        shadedErrorBar(freq_ax,mean(pow(:,kk,:)),std(pow(:,kk,:))./sqrt(size(pow,1)),{'Color',cs(c_ind,:)}),xlim([5 20]), set(gca,'box','off'), set(gca,'XTick',[5 7 10 13 15 20]), ...
            set(gca,'FontName','Arial','FontSize',8), hold on
        c_ind=c_ind+1;
    end
        s_ind = s_ind+1;
end
set(gcf,'PaperSize',[20 20])
suptitle('All Thalamic Power Spectra')
pause(.1)
saveas(gcf,'AllThalFFT.fig')
save('AllThalFFT.pdf')
%%
s_ind=1;     c_ind=1;
%1.6.2 .2.61
cs=[191 206 22; 38 170 73; 237 28 36; 184 48 159; 74 174 230; 50 74 166; 241 78 28; 0 255 0; 128 128 0]./255;
for k=[1 2 3 10 13 14 15 51 52]

    load(strcat('L',num2str(k),'Resultsv5.mat'),'freq','freq_ax','lbl')
    load(strcat('L',num2str(k),'ThChs.mat'))
    CtxChs=1:length(lbl); CtxChs(ThChs)=[]; [~,ll]=min(abs(freq_ax-7)); [~,ul]=min(abs(freq_ax-13));
    pow=freq.fourierspctrm(:,1:length(lbl),:); pow=pow.*conj(pow); 
    [~,fll]=min(abs(freq.freq-2));
    [~,ful]=min(abs(freq.freq-45));
%     mpow=squeeze(mean(pow)); 
%     m=squeeze(mean(mpow,2));
    %pow=pow./repmat(m',[size(pow,1) 1 size(pow,3)]);
    for kk=1:size(pow,2)
        tmp2=pow(:,kk,fll:ful);
        pow(:,kk,:)=pow(:,kk,:)./repmat(mean(tmp2(:)),[size(pow,1) 1 size(pow,3)]);
    end
    %tmpc=pow(:,CtxChs,fll:ful); tmpc=tmpc(:);
    %tmpt=pow(:,ThChs,fll:ful); tmpt=tmpt(:);
    %pow(:,CtxChs,:)=pow(:,CtxChs,:)./sum(tmpc); pow(:,ThChs,:)=pow(:,ThChs,:)./sum(tmpt);
    cpow=squeeze(sum(pow(:,:,ll:ul))); cpow=squeeze(sum(cpow,2));
    tpow=squeeze(sum(pow(:,:,ll:ul))); tpow=squeeze(sum(tpow,2));
    [~,ci]=max(cpow(CtxChs)); ci=CtxChs(ci);
    [~,ti]=max(cpow(ThChs)); ti=ThChs(ti);
        subplot(3,3,s_ind)
        shadedErrorBar(freq_ax,mean(pow(:,ci,:)),std(pow(:,ci,:))./sqrt(size(pow,1)),{'Color','k'}),xlim([4 20]), set(gca,'box','off'), set(gca,'XTick',[4 5 7 10 13 15 20]), ...
            set(gca,'FontName','Arial','FontSize',8), hold on
        shadedErrorBar(freq_ax,mean(pow(:,ti,:)),std(pow(:,ti,:))./sqrt(size(pow,1)),{'Color',cs(c_ind,:)}),xlim([5 20]), set(gca,'box','off'), set(gca,'FontName','Arial','FontSize',8), ...
            set(gca,'XTick',[5 7 10 13 15 20])
        c_ind=c_ind+1;
        s_ind = s_ind+1;
end
set(gcf,'PaperSize',[20 20])
suptitle('Corticothalamic Power Spectra')
pause(.1)
saveas(gcf,'CTFFT.fig')
save('CTFFT.pdf')