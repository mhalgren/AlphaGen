nms={'BW34','Pt1','Pt2','Pt3','Pt5','ChiBi'};

for k=1:length(nms)
    load(strcat(nms{k},'TWResEC.mat'),'agspdma','instfma'); agspdmaec=agspdma; instfmaec=instfma;
    load(strcat(nms{k},'TWResEO.mat'),'agspdma','instfma')
    figure(1)
    subplot(2,2,k)
    hist(agspdma(:,2),100), hold on, hist(agspdmaec(:,2),100)
    figure(2)
    subplot(2,2,k)
    hist(instfma(:,2),100), hold on, hist(instfmaec(:,2),100)
end

    [ky,kx]=ksdensity(agspdma(:,3));
    [~,maxind]=max(ky); spd(k)=kx(maxind);
    [ky,kx]=ksdensity(instfma(:,3));
    [~,maxind]=max(ky); instf(k)=kx(maxind);
    
    
for k=1:length(nms)
    load(strcat(nms{k},'TWResEC.mat'),'agspdma','instfma'); agspdmaec=agspdma; instfmaec=instfma;
    load(strcat(nms{k},'TWResEO.mat'),'agspdma','instfma')

    [ky,kx]=ksdensity(agspdma(:,3)); [kyec,kxec]=ksdensity(agspdmaec(:,3));figure(1), subplot(2,2,k), plot(kx,ky), hold on, plot(kxec,kyec,'r')    
    [pv]=ranksum(agspdma(:,3), agspdmaec(:,3)); title(num2str(pv))
    [ky,kx]=ksdensity(instfma(:,3)); [kyec,kxec]=ksdensity(instfmaec(:,3));figure(2), subplot(2,2,k), plot(kx,ky), hold on, plot(kxec,kyec,'r')
    [pv]=ranksum(instfma(:,3), instfmaec(:,3)); title(num2str(pv))
end

%% Plots frequency histograms of travelling wave speeds and instantaneous frequencies

nms={'Pt1','Pt2','Pt3','Pt5','ChiBi''Pt4'};

for k=1:length(nms)
    if k~=3 && k~=7 && k~=9
    load(strcat(nms{k},'TWResEC.mat'),'agspdma','instfma');
    else
    load(strcat(nms{k},'TWRes.mat'),'agspdma','instfma'); 
    end
    spd=agspdma(:,2);
    if k==6
        spd=spd./2;
    end
    instf=instfma(:,2); spd(spd<0)=[]; instf(instf<0)=[]; [kys,kxs]=ksdensity(spd,0:.01:max(spd)); [kyf,kxf]=ksdensity(instf,0:.01:max(instf));
    figure(1)
    subplot(3,3,k)
    histnorm(spd,0:.01:max(spd)), xlim([0 6])
    title(median(spd)), sp(k)=median(spd);
    %hold on, plot(histXout,pdf(fitdist(spd,'normal'),histXout),'r')
    figure(2)
    subplot(3,3,k)
    histnorm(instf,0:.01:max(instf)), xlim([0 6])
    hold on, plot(kxf,kyf/100,'r'), xlim([5 15])
    title(median(instf)), sf(k)=median(instf);
end
[r,p]=corr(sf',sp','type','spearman');
figure, scatter(sf,sp)
title(strcat('r^2:',{' '},num2str(r.^2),{' '},'p-val:',{' '},num2str(p))), lsline
%% circular plots of TW direction (across time)

nms={'Pt1','Pt2','Pt3','Pt5','ChiBi','Pt4'};
k=5;
for k=1:length(nms)
    if k~=3 && k~=7 && k~=9
    load(strcat(nms{k},'TWResEC.mat'),'ag');
    else
    load(strcat(nms{k},'TWRes.mat'),'ag');
    end
    if k==1
        ag=permute(ag,[2 1 3]); ag=flipdim(ag,1); ag=flipdim(ag,2); %ag=ag(1:5,:,:);
    elseif k==5
        %ag=ag(:,2:end,:);
    elseif k==3
        ag=permute(ag,[2 1 3]);
    elseif k==4
        %ag=ag(:,1:3,:);
    elseif k==6
        ag=flipdim(flipdim(ag,1),2);
    elseif k==7
        %ag=flipdim(ag,2);
    elseif k==8
        ag=flipdim(ag,1);
    end
    tx=gradient(unwrap(squeeze(ag),[],2));
    [~,ty]=gradient(unwrap(squeeze(ag),[],1));
    tx=-tx(:,:,1:(end-1)); ty=ty(:,:,1:(end-1));
    angs=atan2(ty,tx);
    muang=NaN(1,size(angs,3));
    for kk=1:size(angs,3)
        tmp=angs(:,:,kk);
        muang(kk)=circ_mean(tmp(:));
    end
    smu(k)=circ_mean(muang(:));
    figure(1)
    subplot(3,3,k)
    histnorm(muang,20)
    figure(2)
    subplot(3,3,k)
    rose(muang,20)
    figure(3)
    subplot(3,3,k)
    compass(cos(smu(k)),sin(smu(k)))

    %title(circ_median(angs'))
end
figure, compass(cos(smu),sin(smu))
%% circular plots of TW direction (on average)

nms={'Pt1','Pt2','Pt3','Pt5','ChiBi','Pt4'};

for k=1:length(nms)
    if k~=3 && k~=7 && k~=9
    load(strcat(nms{k},'TWResEC.mat'),'agmua');
    else
    load(strcat(nms{k},'TWRes.mat'),'agmua');
    end
    ag=squeeze(agmua(:,:,2));
    if k==1
        ag=permute(ag,[2 1 3]); ag=flipdim(ag,1); ag=flipdim(ag,2); %ag=ag(1:5,:,:);
    elseif k==2||k==5
        %ag=ag(:,2:end,:);
    elseif k==3
        ag=permute(ag,[2 1 3]);
    elseif k==4
        ag=ag(:,1:3,:);
    elseif k==6
        ag=flipdim(flipdim(ag,1),2);
    elseif k==7
        %ag=flipdim(ag,2);
    elseif k==8
        ag=flipdim(ag,1);
    end
    ag(find(ag==0))=NaN;
    tx=-gradient(unwrap(squeeze(ag),[],2));
    [~,ty]=gradient(unwrap(squeeze(ag),[],1));
    angs=atan2(ty,tx);
    figure(1)
    subplot(3,3,k)
    histnorm(angs(:),10)
    figure(2)
    subplot(3,3,k)
    rose(angs(:),10)
    figure(3)
    subplot(3,3,k)
    compass(cos(circ_mean(angs(:))),sin(circ_mean(angs(:))))
    smu(k)=circ_mean(angs(:));
    %title(circ_median(angs'))
end
figure, compass(cos(smu),sin(smu))

%% example (chibi)
    load('ChiBiTWResEC.mat','agmua');
    ag=squeeze(agmua(:,:,2)); ag=flipdim(flipdim(ag,1),2);
    ag(find(ag==0))=NaN;
    tx=-gradient(unwrap(squeeze(ag),[],2));
    [~,ty]=gradient(unwrap(squeeze(ag),[],1));
    angs=atan2(ty,tx); %[x,y]=meshgrid(1:8);
    angmu=atan2(mean(-ty(:)),mean(tx(:)));
    
    figure(1)
    subplot(221)
    imagesc(ag,[-pi/4 pi/4]), axis square, hold on
    quiver(tx,-ty,'LineWidth',4,'Color','k'), hold on, quiver((size(ag,2)./2)+.5,(size(ag,1)./2)+.5,cos(angmu),sin(angmu),'LineWidth',4,'Color','r')%, xlim([2.5 8.5]), ylim([2.5 8.5])
    subplot(222)
    rose(angs(:),10)
    figure(2)
    subplot(2,4,k)
    rose(angs(:),10)
    %title(circ_median(angs'))
%% plot evoked power (of each subject?)
clear
nms={'Pt1','Pt3','Pt5','ChiBi','Pt4'};

figure
for k=1:length(nms)
    load(strcat(nms{k},'TWResEC.mat'),'ecpwrgrd');
    load(strcat(nms{k},'TWResEO.mat'),'eopwrgrd');
    subplot(2,3,k)
    shadedErrorBar(ecpwrgrd.freq,mean(ecpwrgrd.powspctrm),std(ecpwrgrd.powspctrm)./sqrt(size(ecpwrgrd.powspctrm,1)),'b'), hold on
    shadedErrorBar(eopwrgrd.freq,mean(eopwrgrd.powspctrm),std(eopwrgrd.powspctrm)./sqrt(size(eopwrgrd.powspctrm,1)),'r'), xlim([4 25])
    title(size(ecpwrgrd.powspctrm,1))
end

eca=[]; eoa=[];
for k=1:length(nms)
    load(strcat(nms{k},'TWResEC.mat'),'ecpwrgrd');
    load(strcat(nms{k},'TWResEO.mat'),'eopwrgrd');
    [~,fi]=min(abs(ecpwrgrd.freq-4)); [~,fm]=min(abs(ecpwrgrd.freq-25));
    freq_ax=eopwrgrd.freq;
    ecpwrgrd=ecpwrgrd.powspctrm(:,fi:fm);
    ecpwrgrd=ecpwrgrd./repmat(mean(ecpwrgrd,2),[1 size(ecpwrgrd,2)]);
    %ecpwrgrd=ecpwrgrd./mean(ecpwrgrd(:));
    eca=[eca; ecpwrgrd];
    eopwrgrd=eopwrgrd.powspctrm(:,fi:fm);
    eopwrgrd=eopwrgrd./repmat(mean(eopwrgrd,2),[1 size(eopwrgrd,2)]);
    %eopwrgrd=eopwrgrd./mean(eopwrgrd(:));
    eoa=[eoa; eopwrgrd];
end
figure, shadedErrorBar(freq_ax(fi:fm),mean(eca),std(eca)./sqrt(size(eca,1)),'b'), hold on
shadedErrorBar(freq_ax(fi:fm),mean(eoa),std(eoa)./sqrt(size(eoa,1)),'r')

aec=squeeze(mean(eca(:,10:28),2)); aeo=squeeze(mean(eoa(:,10:28),2));
signrank(aec,aeo)

nms={'Pt1','Pt3','Pt5','ChiBi','Pt4'};
eca=[]; eoa=[];
for k=1:length(nms)
    load(strcat(nms{k},'TWResEC.mat'));
    load(strcat(nms{k},'TWResEO.mat'));
    [~,fi]=min(abs(ecpwrgrd.freq-4)); [~,fm]=min(abs(ecpwrgrd.freq-15));
    freq_ax=eopwrgrd.freq;
    if k==1
    eca=[eca; pwrcnv(ecpwrgrd,fi,fm)];
    eca=[eca; pwrcnv(ecpwrstrpr,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrgrd,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrstrpr,fi,fm)];
    elseif k==2
    eca=[eca; pwrcnv(ecpwrgrd,fi,fm)];
    eca=[eca; pwrcnv(ecpwrmo,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrgrd,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrmo,fi,fm)];
    elseif k==3
    eca=[eca; pwrcnv(ecpwrgrd,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrgrd,fi,fm)];
    elseif k==4
    eca=[eca; pwrcnv(ecpwrgrd,fi,fm)];
    eca=[eca; pwrcnv(ecpwrit,fi,fm)];
    eca=[eca; pwrcnv(ecpwrto,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrgrd,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrit,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrto,fi,fm)];
    elseif k==5
    eca=[eca; pwrcnv(ecpwrgrd,fi,fm)];
    eca=[eca; pwrcnv(ecpwrmed,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrgrd,fi,fm)];
    eoa=[eoa; pwrcnv(eopwrmed,fi,fm)];
    end
end
figure, shadedErrorBar(freq_ax(fi:fm),mean(eca),std(eca)./sqrt(size(eca,1)),'b'), hold on
shadedErrorBar(freq_ax(fi:fm),mean(eoa),std(eoa)./sqrt(size(eoa,1)),'r'), xlim([4 15])

signrank(mean(eca(:,10:28),2),mean(eoa(:,10:28),2))