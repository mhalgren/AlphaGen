%% FIND ACTUAL INTERCONTACT DISTANCE!!!
clear
wdir='';
addpath(genpath(strcat(wdir,'fieldtrip-20150121')))
cd(strcat(wdir,'AlphaLaminar'))
addpath(genpath('Scripts'))
addpath(genpath('PGData'))
cd('TW')
addpath(genpath('circ_stat'))
addpath(genpath('Monkeys'))
clear
for ec=[1 0]
    %load data
    if ec==1
        load('ChibiSess1Data.mat')
    else
        load('ChibiSess1DataEO.mat')
    end
origdata=data; clear data

frq=[7 13];

    if ec==1
        load('ChiBiArt.mat')
    else
        load('ChiBiArtEO.mat')
    end


for f=1:size(frq,1)
%% medial strip analysis
cfg=[];
cfg.channel=origdata.label([121:128]);
data=ft_preprocessing(cfg,origdata);

trl=[];

if ec==1
    ecpwrmed=evokedpwr(data,trl,tmp);
else
    eopwrmed=evokedpwr(data,trl,tmp);
end

cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[frq(f,1) frq(f,2)];

adata=ft_preprocessing(cfg,data);
cfg=[]; cfg.hilbert='abs';
adatap=ft_preprocessing(cfg,adata);
cfg=[];
cfg.hilbert='angle';
adata=ft_preprocessing(cfg,adata);

cfg=[];
cfg.length=3;
adata=ft_redefinetrial(cfg,adata);
adatap=ft_redefinetrial(cfg,adatap);

art_smpl=tmp.artfctdef.visual.artifact;
bad_trials=zeros(1,length(art_smpl));
for kk=1:size(art_smpl,1);
    smp=round((art_smpl(kk,1)+art_smpl(kk,2))./2);
    for tind=1:length(adata.sampleinfo)
        if smp<=adata.sampleinfo(tind,2) && smp>=adata.sampleinfo(tind,1)
        bad_trials(kk)=tind;
        break
        end
    end
end

bad_trials(find(bad_trials==0))=[];

cfg=[];
trs=1:length(adata.trial);
trs(bad_trials)=[];
cfg.trials=trs;
adata=ft_redefinetrial(cfg,adata);
adatap=ft_redefinetrial(cfg,adatap);

%% order into matrix of chanxchanxtime of phs diff from mean
a=fieldtrip2mat(adata); ap=fieldtrip2mat(adatap);
adiff=circ_dist(a,repmat(circ_mean(a),[size(a,1) 1]));
p=nanmean(ap);
[~,ind]=sort(p);
ind=ind(round(length(ind)*.2):end);
ind=sort(ind); ind(1)=[]; ind(end)=[];

adiff=adiff(:,ind);
a=a(:,ind);
med=a;

%% grid analysis
cfg=[];
cfg.channel={'75','76','77','78','79','87','88','89','90','91','100','101','102','103','104','110','111','112','113','114'};
data=ft_preprocessing(cfg,origdata);

sdata=plateaufilthilb(data,origdata.fsample,frq(f,:),.15);

if ec==1
    ecpwrgrd=evokedpwr(data,trl,tmp);
else
    eopwrgrd=evokedpwr(data,trl,tmp);
end
%% filter into bands

cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[frq(f,1) frq(f,2)];

adata=ft_preprocessing(cfg,data);
cfg=[]; cfg.hilbert='abs';
adatap=ft_preprocessing(cfg,adata);
cfg=[];
cfg.hilbert='angle';
adata=ft_preprocessing(cfg,adata);

cfg=[];
cfg.length=3;
adata=ft_redefinetrial(cfg,adata);
adatap=ft_redefinetrial(cfg,adatap);
sdata=ft_redefinetrial(cfg,sdata);

art_smpl=tmp.artfctdef.visual.artifact;
bad_trials=zeros(1,length(art_smpl));
for kk=1:size(art_smpl,1);
    smp=round((art_smpl(kk,1)+art_smpl(kk,2))./2);
    for tind=1:length(adata.sampleinfo)
        if smp<=adata.sampleinfo(tind,2) && smp>=adata.sampleinfo(tind,1)
        bad_trials(kk)=tind;
        break
        end
    end
end

bad_trials(find(bad_trials==0))=[];

cfg=[];
trs=1:length(adata.trial);
trs(bad_trials)=[];
cfg.trials=trs;
adata=ft_redefinetrial(cfg,adata);
adatap=ft_redefinetrial(cfg,adatap);
sdata=ft_redefinetrial(cfg,sdata);
%% order into matrix of chanxchanxtime of phs diff from mean
a=fieldtrip2mat(adata); ap=fieldtrip2mat(adatap);
bchs=[];
gchs=1:size(a,1); gchs(bchs)=[];
adiff=circ_dist(a,repmat(circ_mean(a(gchs,:)),[size(a,1) 1]));
adiff(bchs,:)=NaN(length(bchs),size(adiff,2));
a(bchs,:)=NaN(length(bchs),size(a,2));
ap(bchs,:)=NaN(length(bchs),size(ap,2));
%adiff=circ_dist(a,repmat(a(20,:),[size(a,1) 1]));
%adiff=circ_dist(a,repmat(circ_mean(a([1 9 17 25 2 10 18 26 3 11 19 27],:)),[size(a,1) 1]));
%adiff=circ_dist(a,repmat(circ_mean(a([17 25 18 26 19 27],:)),[size(a,1) 1]));

agdiff=zeros(4,5,size(adiff,2));
ag=agdiff; 
if f==3
    apord=agdiff;
end
instf=agdiff(:,:,1:(size(adiff,2)-1));
instf_unord=ft_preproc_medianfilter(diff(unwrap(fieldtrip2mat(sdata),[],2),[],2),10);

a_ind=1;
for k=1:4
    for kk=1:5
        disp(kk)
        agdiff(k,kk,:)=adiff(a_ind,:);
        ag(k,kk,:)=a(a_ind,:);
        instf(k,kk,:)=instf_unord(a_ind,:);
        if f==2
        apord(k,kk,:)=ap(a_ind,:);
        end
        a_ind=a_ind+1;
    end
end

tx=gradient(unwrap(squeeze(ag),[],2));
[~,ty]=gradient(unwrap(squeeze(ag),[],1));
tx=tx(:,:,1:(end-1)); ty=ty(:,:,1:(end-1));
angs=atan2(ty,tx);
agspd=instf./sqrt(tx.^2+ty.^2);

%% find top 20% of alpha epochs
p=nanmean(ap);
[~,ind]=sort(p);
ind=ind(round(length(ind)*.8):end);
ind=sort(ind); ind(1)=[]; ind(end)=[];

agdiff=agdiff(:,:,ind);
agspd=agspd(:,:,ind);
angs=angs(:,:,ind);
instf=instf(:,:,ind);
ag=ag(:,:,ind);

agmu=circ_mean(agdiff,[],3);

for k=1:size(agspd,3)
    tmp1=agspd(:,:,k); tmp2=instf(:,:,k);
    tmpang=angs(:,:,k);tmpang=tmpang(:); tmpang=tmpang(~isnan(tmpang));    
    instfm(k)=nanmedian(tmp2(:));
    agspdm(k)=nanmedian(tmp1(:));
    instangs(k)=circ_mean(tmpang(:));
end

agspdm=agspdm.*origdata.fsample./100;
instfm=instfm.*origdata.fsample./(2*pi);

agspdma(:,f)=agspdm;
instfma(:,f)=instfm;

for k=1:size(agspd,1)
    for kk=1:size(agspd,2)
        [tmpr,tmpp]=corrcoef(squeeze(agspd(k,kk,:)),squeeze(instf(k,kk,:)));
        r(k,kk,f)=tmpr(1,2);
        pv(k,kk,f)=tmpp(1,2);
    end
end
%imagesc(agmu)
% imagesc(circ_mean(agdiff,[],3))
% 
% x=unwrap(circ_mean(agdiff,[],3),[],2);
% x=zscore(x,[],2);
% imagesc(x)
% imagesc()

%% save in overall matrices
agmua(:,:,f)=agmu;
meda(:,:,f)=med;
if f==1
    agdiffd=agdiff;
end
end
if ec
save('ChiBiTWResEC.mat','agspdma','instfma','r','pv','agmua','meda','ag','apord','agdiff','agdiffd','ecpwrmed','ecpwrgrd','-v7.3')
else
    save('ChiBiTWResEO.mat','agspdma','instfma','r','pv','agmua','meda','ag','apord','agdiff','agdiffd','eopwrmed','eopwrgrd','-v7.3')
end
clear
end