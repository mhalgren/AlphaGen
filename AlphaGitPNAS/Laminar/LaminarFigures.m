snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r'};
clear
figure
for k=1:3
    load(strcat('LoopOut10/',snm{k},'_AlphaFreqNu.mat'))
%     if k~=3
%     x=frqalph.fourierspctrm.*conj(frqalph.fourierspctrm);
%     else
        x=frqalph.powspctrm;
   % end
    [~,xmuchan]=max(mean(squeeze(mean(x(:,:,find(frqalph.freq==7):find(frqalph.freq==13)),3)))); x=squeeze(x(:,xmuchan,:));
    x=x./repmat(mean(x(:)),[size(x,1) size(x,2)]);
    pow(k,:) = squeeze(mean(x)); sem(k,:) = squeeze(std(x))./size(x,1);
    shadedErrorBar(frqalph.freq,squeeze(mean(x)), squeeze(std(x))./size(x,1),'lineprops',cs{k}), xlim([4 25]), hold on
end
%
fq = frqalph.freq;
%
title('Normalized CSD Power Spectra'), ylabel('Power'), xlabel('Frequency (Hz)'),  set(gca,'box','off','FontName','Arial','FontSize',8)
%%
snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r','k'}; sys=[20 10 10];
cl=[-409.6 409.6;-80 80; -1.4 1.4];
figure
for k=1:3
    load(strcat('LoopOut10/',snm{k},'_AvgAlpha.mat'))
    subplot(2,2,k)
    x=alpha_avg.avg; 
%     
%     gx=-2:2;
%     gaus2d=zeros(length(gx));
%     sx=.75;
%     sy=sys(k)%10^-100;
%     
%     for xi=1:length(gx)
%         for yi=1:length(gx)
%             gaus2d(xi,yi)=exp(-((gx(xi)^2)/(2*sx^2)+(gx(yi)^2)/(2*sy^2)));
%         end
%     end
%     x=conv2(x,gaus2d,'same'); x=x./sum(gaus2d(:));
    
    colormap jet

    imagesc(linspace(-250,250,length(x)),1:size(x,1),x), set(gca,'box','off')
    avgcsd{k} = x;
    tx_csd{k} = linspace(-250,250,length(x));
    switch k
        case 1
        set(gca,'YTick',[2 5.5 10.5 14.5 18.5 22]-.5,'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
        case 2
        set(gca,'YTick',[2 4.5 8.5 11.5 13 16]-.5,'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
        case 3
        set(gca,'YTick',[1 3.5 9 14.5 18 21]-1.5,'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
    end
    caxis(cl(k,:))
end
title('Alpha CSD Profiles'), ylabel('Estimated Cortical Layer'), xlabel('Time (ms)')

snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r','k'}; sys=[20 10 10];
figure
for k=1:2
    load(strcat('LoopOut10/',snm{k},'_AvgMUA.mat')), load(strcat('LoopOut10/',snm{k},'_AlphaMUABl.mat'))
    subplot(2,2,k)
    if k==2
        mua_avg.avg=mua_avg.avg(1:17,:);
        muamu=muamu(1:17);
    end
    x=mua_avg.avg-repmat(muamu,[1 size(mua_avg.avg,2)]); 
    
    gx=-2:2;
    gaus2d=zeros(length(gx));
    sx=.75;
    sy=sys(k);
    
    for xi=1:length(gx)
        for yi=1:length(gx)
            gaus2d(xi,yi)=exp(-((gx(xi)^2)/(2*sx^2)+(gx(yi)^2)/(2*sy^2)));
        end
    end
    x=conv2(x,gaus2d,'same'); x=x./sum(gaus2d(:));
    imagesc(linspace(-250,250,length(x)),1:size(x,1),x), set(gca,'box','off')
    avgmua{k} = x;
    tx_mua{k} = linspace(-250,250,length(x));
    switch k
        case 1
        set(gca,'YTick',[2 5.5 10.5 14.5 18.5 22]-.5,'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
        case 2
        set(gca,'YTick',[2 4.5 8.5 11.5 13 16]-.5,'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
        case 3
        otherwise
    end
    caxis([-.4 .4])
end

title('Alpha MUA Profiles'), ylabel('Estimated Cortical Layer'), xlabel('Time (ms)')

%%
snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r','k'}; sys=[20 10 10];
figure
for k=1:2
    load(strcat('LoopOut10/',snm{k},'_AvgAlpha.mat')), load(strcat('LoopOut10/',snm{k},'_AlphaModInd.mat'),'hmu')
    subplot(2,2,k)
    x=hgp_avg.avg-repmat(hmu,[1 size(hgp_avg.avg,2)]); 
%     
%     gx=-2:2;
%     gaus2d=zeros(length(gx));
%     sx=.75;
%     sy=sys(k);
%     
%     for xi=1:length(gx)
%         for yi=1:length(gx)
%             gaus2d(xi,yi)=exp(-((gx(xi)^2)/(2*sx^2)+(gx(yi)^2)/(2*sy^2)));
%         end
%     end
%     x=conv2(x,gaus2d,'same'); x=x./sum(gaus2d(:));
    imagesc(linspace(-250,250,length(x)),1:size(x,1),x), set(gca,'box','off')
    hgp_avgprof{k} = x; hgp_tx{k} = linspace(-250,250,length(x));
    switch k
        case 1
        set(gca,'YTick',[1 4 9 14.5 16 21],'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
        case 2
        set(gca,'YTick',[2 4.5 8.5 11.5 13 16]-.5,'YTickLabel',{'I','II','III','IV','V','VI'},'XTick',[-250 -125 0 125 250])
        case 3
        otherwise
    end
    %caxis([-.4 .4])
end
title('Alpha HGP Profiles'), ylabel('Estimated Cortical Layer'), xlabel('Time (ms)')
% %%
% figure
% snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r'};
% for k=1:2
%     load(strcat('LoopOut9/',snm{k},'_AvgAlpha.mat')),load(strcat('LoopOut9/',snm{k},'_AlphaModIndMUA.mat'))
%     z = sqrt(2) * erfcinv((.05/23)*2);
%     [~,am]=max(mean(alpha_avg.avg.^2,2));
%     plot((squeeze(ModIndZma(am,:))),cs{k}), hold on, plot(1:23,ones(1,23).*z,'r'), set(gca,'XTick',[1 4 9 14 16 21],'XTickLabel',{'I','II','III','IV','V','VI'})
% end
% title('Modulation Index'), ylabel('Z-Score'), xlabel('Estimated Cortical Layer'),  set(gca,'box','off'), ylim([-1 30])

figure
colormap jet
snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r'};
for k=1:2
    load(strcat('LoopOut10/',snm{k},'_AlphaModIndMUA.mat'))
    if k==2
        ModIndZma=ModIndZma(:, 1:17);
    end
    nch=size(ModIndZma,1);
    MI=ModIndZma(1:nch,1:end); sig=MI>(sqrt(2)*erfcinv((.05/length(MI(:)))*2));
    hgp_mi{k} = MI';
    subplot(1,2,k)
    contourf(MI',40,'linecolor','none'), set(gca,'box','off'), axis equal, hold on
    contour(sig',1,'linecolor','k','linew',6)
    
    switch k
        case 1
        set(gca,'XTick',[2 4.5 9.5 14 18 22],'XTickLabel',{'I','II','III','IV','V','VI'},'YTick',[1 4.5 9.5 14 18 22],'YTickLabel',{'I','II','III','IV','V','VI'})
        case 2
        set(gca,'XTick',[2 4.5 8.5 11.5 13 16]-1,'XTickLabel',{'I','II','III','IV','V','VI'},'YTick',[2 4.5 8.5 11.5 13 16],'YTickLabel',{'I','II','III','IV','V','VI'})
        case 3
        otherwise
    end
    
    if k==1
        caxis([0 15])
    else
        caxis([0 15])
    end
end
title('Modulation Index'), ylabel('Estimated Cortical Layer'), xlabel('Estimated Cortical Layer')
%%
hgp_mi = {};
figure
colormap jet
snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r'};
for k=1:2
    load(strcat('LoopOut10/',snm{k},'_AlphaModInd.mat'))
    if k==2
        ModIndZa=ModIndZa(:, 1:17);
    end
    nch=size(ModIndZa,1);
    MI=ModIndZa(1:nch,1:end); sig=MI>(sqrt(2)*erfcinv((.05/length(MI(:)))*2));
    hgp_mi{k} = MI';
    subplot(1,2,k)
    contourf(MI',40,'linecolor','none'), set(gca,'box','off'), axis equal, hold on
    contour(sig',1,'linecolor','k','linew',6)
    
    switch k
        case 1
        set(gca,'XTick',[2 4.5 9.5 14 18 22],'XTickLabel',{'I','II','III','IV','V','VI'},'YTick',[1 4.5 9.5 14 18 22],'YTickLabel',{'I','II','III','IV','V','VI'})
        case 2
        set(gca,'XTick',[2 4.5 8.5 11.5 13 16]-1,'XTickLabel',{'I','II','III','IV','V','VI'},'YTick',[2 4.5 8.5 11.5 13 16],'YTickLabel',{'I','II','III','IV','V','VI'})
        case 3
        otherwise
    end
    
    if k==1
        caxis([0 25])
    else
        caxis([0 25])
    end
end
title('HGP Modulation Index'), ylabel('Estimated Cortical Layer'), xlabel('Estimated Cortical Layer')
%%
figure
snm={'Pt2','Pt3','Pt1'}; cs={'b','m','r'}; cf=[3000 300]; mm=[11 9]; sf=[20 10];
for k=1:2
    load(strcat('LoopOut10/',snm{k},'_AvgAlpha.mat'))
    if k~=3
        load(strcat('LoopOut10/',snm{k},'_AvgMUA.mat'));
    end
    subplot(1,2,k)
    [~,am]=max(mean(alpha_avg.avg.^2,2));
    shadedErrorBar(alpha_avg.time,(squeeze(alpha_avg.avg(am,:))),(squeeze(sqrt(alpha_avg.var(am,:))./sqrt(length(alpha_avg.cfg.previous.trl)))),'lineprops','g'), ...
        ylabel('Power (a.u.)'), xlabel('Time (s)'),  set(gca,'box','off'), xlim([-.25 .25]), set(gca,'ytick',[],'xtick',[-.25 0 .25]), hold on
    if k~=3
        shadedErrorBar(mua_avg.time,smooth((squeeze(mua_avg.avg(mm(k),:))-mean(mua_avg.avg(mm(k),:))).*cf(k),sf(k)),(squeeze(sqrt(cf(k).*mua_avg.var(mm(k),:))./sqrt(length(mua_avg.cfg.previous.trl)))),'lineprops','m'), hold on
    end
end
suptitle('CSD vs. MUA')

%% plot coherence
snm={'Pt2','Pt3'}; cs={'b','m','r','k'};
figure
for k=1:2
    load(strcat('LoopOut6/',snm{k},'_AlphaCoh.mat'))
    subplot(1,2,k)
    cohshff=squeeze(mean(cohall(:,:,11:23,:),3));
    cohshffmu=mean(cohshff,3); cohshffstd=std(cohshff,[],3);
    acoh=squeeze(mean(coh.cohspctrm(:,:,11:23),3));
    acohz=(acoh-cohshffmu)./cohshffstd;
    acohz=acohz(1:23,49:end); asig=acohz>3.7332;
     acoh=acoh(1:23,49:end);
    contourf(acoh,40,'linecolor','none'), set(gca,'box','off'), set(gca,'XTick',[1 4 9 14 16 21],'XTickLabel',{'I','II','III','IV','V','VI'},'YTick',[1 4 9 14 16 21],'YTickLabel',{'I','II','III','IV','V','VI'}), axis equal
    hold on
    contour(asig,1,'linecolor','k','linew',2)
    caxis([0 .35])
end
title('Alpha CSD Profiles'), ylabel('Estimated Cortical Layer'), xlabel('Time (ms)'),  set(gca,'box','off')
%%

ad = ([2 5.5 10.5 14.5 18.5 22]) + ([2 4.5 8.5 11.5 13 16]) + ([1 3.5 9 14.5 18 21]); ad=ad./3;

figure
load('LoopOut10/Pt2_AlphaFreqNu.mat')
nch=length(frqalph.label);
pow=frqalph.powspctrm; pow=pow./repmat(mean(pow(:)),[size(pow,1) size(pow,2) size(pow,3)]);
powcsd{1}=squeeze(mean(mean(pow(:,1:nch,31:61),3)));
semcsd{1}=squeeze(std(mean(pow(:,1:nch,31:61),3)))./sqrt(size(frqalph.powspctrm,1));
shadedErrorBar(1:nch,squeeze(mean(mean(pow(:,1:nch,31:61),3))),squeeze(std(mean(pow(:,1:nch,31:61),3)))./sqrt(size(frqalph.powspctrm,1)),'lineprops','b'), xlim([1 24]), set(gca,'box','off')
ylabel('Normalized Power'), xlabel('Estimated Cortical Layer')
hold on

load('LoopOut10/Pt3_AlphaFreqNu.mat') 
nch=length(frqalph.label);
pow=frqalph.powspctrm; pow=pow./repmat(mean(pow(:)),[size(pow,1) size(pow,2) size(pow,3)]);
powcsd{2}=squeeze(mean(mean(pow(:,1:nch,31:61),3)));
semcsd{2}=squeeze(std(mean(pow(:,1:nch,31:61),3)))./sqrt(size(frqalph.powspctrm,1));
shadedErrorBar(1:nch,squeeze(mean(mean(pow(:,1:nch,31:61),3))),squeeze(std(mean(pow(:,1:nch,31:61),3)))./sqrt(size(frqalph.powspctrm,1)),'lineprops','m'), xlim([1 24]), set(gca,'box','off')
hold on
load('LoopOut10/Pt1_AlphaFreqNu.mat')
nch=length(frqalph.label);
pow=frqalph.powspctrm; pow=pow./repmat(mean(pow(:)),[size(pow,1) size(pow,2) size(pow,3)]);
powcsd{3}=squeeze(mean(mean(pow(:,1:nch,31:61),3)));
semcsd{3}=squeeze(std(mean(pow(:,1:nch,31:61),3)))./sqrt(size(frqalph.powspctrm,1));
shadedErrorBar([1:nch].*(1.75/1.5),squeeze(mean(mean(pow(:,1:nch,31:61),3))),squeeze(std(mean(pow(:,1:nch,31:61),3)))./sqrt(size(frqalph.powspctrm,1)),'lineprops','r'), xlim([1 24]), set(gca,'box','off'), set(gca,'XTick',ad,'XTickLabel',{'I','II','III','IV','V','VI'})
hold on

title('Alpha Power vs. Cortical Depth')

