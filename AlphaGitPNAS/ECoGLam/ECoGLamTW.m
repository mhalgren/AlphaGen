%%
clear
load('Pt2_EcogLamData.mat')
load('Pt2_ELTW_Arts.mat')

bad_trials=[82 108:174];

comb_dat.trial{1}(14,:)=(2./3).*comb_dat.trial{1}(13,:) + (1./3).*comb_dat.trial{1}(16,:);
comb_dat.trial{1}(15,:)=(1./3).*comb_dat.trial{1}(13,:) + (2./3).*comb_dat.trial{1}(16,:);

origdata=comb_dat;

data=comb_dat; clear comb_dat

cfg=[];
cfg.lpfilter='yes';
cfg.lpfreq=30;
data=ft_preprocessing(cfg,data);
cfg=[];
cfg.bsfilter='yes';
cfg.bsfreq=[48 52];
data=ft_preprocessing(cfg,data);
cfg=[];
cfg.hpfilter='yes';
cfg.hpfreq=2;
data=ft_preprocessing(cfg,data);

cfg=[];
cfg.channel=data.label(25:64);
data=ft_preprocessing(cfg,data);

sdata=plateaufilthilb(data,data.fsample,[7 13],.15);

cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[7 13];
b=ft_preprocessing(cfg,data);

cfg=[];
cfg.hilbert='abs';
bp=ft_preprocessing(cfg,b);

cfg=[];
cfg.hilbert='angle';
b=ft_preprocessing(cfg,b);

cfg=[]; cfg.length=3; 
b=ft_redefinetrial(cfg,b);
bp=ft_redefinetrial(cfg,bp);
sdata=ft_redefinetrial(cfg,sdata);

cfg=[]; cfg.trials=1:length(b.trial); cfg.trials(bad_trials)=[];
b=ft_redefinetrial(cfg,b);
bp=ft_redefinetrial(cfg,bp);
sdata=ft_redefinetrial(cfg,sdata);

b=fieldtrip2mat(b);
bp=mean(fieldtrip2mat(bp));
[~,ind]=sort(bp);
ind=ind(round(length(ind)*.8):end);
instf_unord=ft_preproc_medianfilter(diff(unwrap(fieldtrip2mat(sdata),[],2),[],2),10);

%%
cfg=[];
cfg.channel=origdata.label(1:18);
ldata=ft_preprocessing(cfg,origdata);

for k=1:length(ldata.trial)
c=ldata.trial{k}; c=pg2csdv3(c,.15); ldata.trial{k}(1:16,:)=c;
end

cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[8 12];
lb=ft_preprocessing(cfg,ldata);

cfg=[];
cfg.hilbert='angle';
lb=ft_preprocessing(cfg,lb);

cfg=[]; cfg.length=3; lb=ft_redefinetrial(cfg,lb);

cfg=[]; cfg.trials=1:length(lb.trial); cfg.trials(bad_trials)=[];
lb=ft_redefinetrial(cfg,lb);

lb=fieldtrip2mat(lb);

%%

%use if thresholding
b=b(:,ind);
lb=lb(:,ind);
instf_unord=instf_unord(:,ind);

bd=circ_dist(b,repmat(circ_mean(b),[size(b,1) 1]));

lbd=circ_dist(lb,repmat(circ_mean(b),[size(lb,1) 1]));
lbt=circ_mean(lbd,[],2);

b_ord=NaN(8,5,length(b));
b_raw=NaN(8,5,length(b));
instf=b_ord;

for k=1:40
    i=str2num(data.label{k}(3:4));
    disp(i)
    [x,y]=ind2sub([8,5],i);
        disp([x,y])
    b_ord(x,y,:)=bd(k,:);
    b_raw(x,y,:)=b(k,:);
    instf(x,y,:)=instf_unord(k,:);
end

b_ord=flipud(b_ord); b_raw=flipud(b_raw); instf=flipud(instf);

bt=circ_mean(b_ord,[],3);

%%
tx=-gradient(unwrap(squeeze(b_raw),[],2));
[~,ty]=gradient(unwrap(squeeze(b_raw),[],1));
b_raw=b_raw(:,:,1:(end-1)); %b_raw=b_raw(:,:,1:(end-1));
angs=atan2(ty,tx);
agspd=instf./sqrt(tx.^2+ty.^2);

for k=1:size(angs,3)
    tmp1=agspd(:,:,k); tmp2=instf(:,:,k);
    tmpang=angs(:,:,k);tmpang=tmpang(:); tmpang=tmpang(~isnan(tmpang));
    instfm(k)=nanmedian(tmp2(:));
    agspdm(k)=nanmedian(tmp1(:));
    instangs(k)=circ_mean(tmpang(:));
    pgd(k)=circ_r(tmpang(:));
    disp(k./size(angs,3))
%     tmpang=angs(:,:,k);tmpang=tmpang(:); tmpang=tmpang(~isnan(tmpang));
%     instangs(k)=circ_mean(tmpang(:)); 
%     instcoh(k)=circ_r(tmpang(:));
end

agspdm=agspdm.*data.fsample./100;
instfm=instfm.*data.fsample./(2*pi);
%imagesc(bt,[-pi/5 pi/5]), colormap jet

figure(1)
histnorm(agspdm,100), box off, set(gca,'YTick',[0 1 2 3],'XTick',[.2 1 2]), set(gcf,'PaperSize',[11 6.5])

figure(2)
[tout, rout] = rose(instangs,50);
polar(tout, rout);
[xout, yout] = pol2cart(tout, rout);
set(gca, 'nextplot', 'add');
fill(xout, yout, 'b');
set(gcf,'PaperSize',[6 6])

rticks(104449*[.01:.01:.06])

btt=flipud(bt);

%load('Pt2Coord.mat')
load('Pt2ElecCoordv2.mat')

lcoor=ones(21,2).*832;
lcoor(:,2)=linspace(165,600,21);

figure(3)
ti=imread('Pt2Brn2.jpg');
imagesc(ti), hold on
scatter(elec_coord(:,1),elec_coord(:,2),250,btt(:),'filled'), caxis([-pi/4 pi/4]), ...
    axis equal, box off, axis off, hold on, scatter(380,338,250,lbt(1),'filled','d'), hold on%scatter(380,338,150,'k','filled','d'), hold on
%scatter(lcoor(:,1),lcoor(:,2),250,lbt(:),'filled'), caxis([-pi/4 pi/4]), ...
    %axis equal, box off, axis off, hold on, scatter(380,338,150,'k','filled','d'), hold on
colormap jet
set(gcf,'PaperSize',[6 6])


%%
ti=imread('Pt2_reco.jpg');
subplot(2,1,1)
imagesc(ti), hold on
scatter(coord(:,1),825-coord(:,2),250,btt(:),'filled'), caxis([-pi/4.5 pi/4.5]), ...
    axis equal, box off, axis off, hold on, scatter(560,360,150,'k','filled','d'), colormap jet

ti=imread('Pt2_reco.jpg');
subplot(2,1,1)
imagesc(ti), hold on
scatter(coord(:,1),825-coord(:,2),250,btt(:),'filled'), caxis([-pi/4.5 pi/4.5]), axis equal, box off, axis off, hold on, scatter(560,360,150,'k','filled','d'), colormap jet
subplot(2,1,2)
imagesc(lbt,[-pi/4.5 pi/4.5])

figure
for k=513:size(b_raw,3)
%     subplot(2,1,1)
    imagesc(squeeze(b_raw(:,:,k)),[-pi pi])
    pause(.5)
        title(k)
        drawnow
        %pause(.5)
%     subplot(2,1,2)
%     imagesc(c(:,(k-512):(k+512)))
%     drawnow
end