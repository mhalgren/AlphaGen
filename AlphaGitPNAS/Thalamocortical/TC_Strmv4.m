function [alpha_td,lbl,freq,coh,cohz,cohreal,cohall,pdc,psi,freq_ax,r,pv,alpha_var,ntrs,mi,mod_ind,SR] = TC_Strmv4(data,btrs,ThChs)

lbl=data.label;
freq=[];
coh=[];
pdc=[];
psi=[];
freq_ax=[];
r=[];
pv=[];
%% reject machine noise
cfg=[];
cfg.bsfilter='yes';
bsmax=(data.fsample./2)-rem((data.fsample./2),50);
bsmax=50:50:bsmax; bsmax=bsmax'; bsmax=[bsmax-2 bsmax+2];
cfg.bsfreq=bsmax;%CHANGE
cfg.bsinstabilityfix='reduce';
data=ft_preprocessing(cfg,data);
% cfg=[];
% cfg.bpfilter='yes';
% cfg.bpfreq=[250 500];%CHANGE
% hhgp=ft_preprocessing(cfg,data);
% cfg=[];
% cfg.hilbert='abs';
% hhgp=ft_preprocessing(cfg,hhgp);
% for k=1:length(hhgp.label)
% hhgp.label{k}=['h' hhgp.label{k}];
% hhgp.trial{1}(k,:)=smooth(hhgp.trial{1}(k,:),20);
% end
if data.fsample>=129
cfg=[];
cfg.bpfilter='yes';
if data.fsample<=257
    cfg.bpfreq=[70 120];
else
cfg.bpfreq=[70 190];%CHANGE
end
lhgp=ft_preprocessing(cfg,data);
cfg=[];
cfg.hilbert='abs';
lhgp=ft_preprocessing(cfg,lhgp);
for k=1:length(lhgp.label)
lhgp.label{k}=['l' lhgp.label{k}];
lhgp.trial{1}(k,:)=smooth(lhgp.trial{1}(k,:),20);
end
end
cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[7 13];
cfg.bpinstabilityfix='reduce';
alph=ft_preprocessing(cfg,data);
cfg=[];
cfg.hilbert='abs';
alph_amp=ft_preprocessing(cfg,alph);
ap=alph_amp.trial{1};
[xca,lga]=xcorr(ap',5.*data.fsample);
cfg=[];
cfg.hilbert='angle';
alph_phs=ft_preprocessing(cfg,alph);
cfg=[];
cfg.hpfilter='yes';
cfg.hpfreq=2;
cfg.hpinstabilityfix='reduce';
data=ft_preprocessing(cfg,data);
cfg=[];
%GIVE LHGP AND HHGP NEW LABELS!!!
if data.fsample>=129
data=ft_appenddata(cfg,data,lhgp);
end
%data=ft_appenddata(cfg,data,lhgp,hhgp);
%% reject trials w/ noise
cfg=[];
cfg.length=4;
cfg.overlap=0;
data=ft_redefinetrial(cfg,data);
origalph=ft_redefinetrial(cfg,alph);
alph_amp=ft_redefinetrial(cfg,alph_amp);
if data.fsample>=129
lhgp=ft_redefinetrial(cfg,lhgp);
alph_phs=ft_redefinetrial(cfg,alph_phs);
end
cfg=[];
trs=1:length(data.trial);
trs(btrs)=[];
cfg.trials=trs;
data=ft_redefinetrial(cfg,data);
origalph=ft_redefinetrial(cfg,origalph);
alph_amp=ft_redefinetrial(cfg,alph_amp);
if data.fsample>=129
lhgp=ft_redefinetrial(cfg,lhgp);
alph_phs=ft_redefinetrial(cfg,alph_phs);
end
%% split into epochs
cfg=[];
cfg.length=2;
cfg.overlap=0;
data=ft_redefinetrial(cfg,data);
origalph=ft_redefinetrial(cfg,origalph);
alph_amp=ft_redefinetrial(cfg,alph_amp);
if data.fsample>=129
lhgp=ft_redefinetrial(cfg,lhgp);
alph_phs=ft_redefinetrial(cfg,alph_phs);
end
%%
orig_data=data;
a=squeeze(mean(squeeze(fieldtrip2mat_epochs(alph_amp)),2));

%% only use alpha trials
data=orig_data;

% HERE IS WHERE I DECIDE WHICH CHANNEL'S ALPHA TO USE FOR PICKING HIGH
% ALPHA- ALL CHANS OR JUST THALAMIC???

atmp=squeeze(mean(a(ThChs,:)));%atmp=squeeze(mean(a));
[~,ind]=sort(atmp); ind=ind(round(length(ind).*.8):end);

%% reject non alpha epochs
cfg=[];
cfg.trials=ind;
data=ft_redefinetrial(cfg,data);
alph=ft_redefinetrial(cfg,origalph);
if data.fsample>=129
lhgp=ft_redefinetrial(cfg,lhgp);
alph_phs=ft_redefinetrial(cfg,alph_phs);
end
%% calculate sharpness ratio
SR=SharpnessRatio(origalph);
%% compute Tort's MI
if data.fsample>=129
    nChans=length(lhgp.label);
mi=zeros(nChans,nChans,36);
mod_ind=zeros(nChans,nChans);
nPerms=200;
for chan = 1:nChans
disp(strcat('Now processing channel',{' '},num2str(chan),{' '},'out of',{' '},num2str(nChans)))

%calculate MI
allphase=fieldtrip2mat(alph_phs);
allamplitude=fieldtrip2mat(lhgp);
    phase=squeeze(allphase(chan,:));
    real_pac=average_envelope_versus_phase_v3(allamplitude,phase);
    tmp=zeros(size(allamplitude,1),36,nPerms);
    for nperm=1:nPerms
        shuffpoint=randi(round(length(phase).*.6))+round(length(phase).*.2);
        shuffphase=squeeze(phase([shuffpoint:end 1:(shuffpoint-1)]));
        tmp(:,:,nperm)=average_envelope_versus_phase_v3(allamplitude,shuffphase);
        %normalize tmp
    end
    mi(chan,:,:)=squeeze((real_pac-mean(tmp,3))./std(tmp,[],3));
    tmp=tmp./repmat(sum(tmp,2),[1 size(tmp,2) 1]);
    real_pac=real_pac./repmat(sum(real_pac,2),[1 size(real_pac,2)]);
    real_pac2=squeeze((log(36)+sum(real_pac.*log(real_pac),2))./log(36));
    tmp2=squeeze((log(36)+sum(tmp.*log(tmp),2))./log(36));
    mod_ind(chan,:)=squeeze((real_pac2-mean(tmp2,2))./std(tmp2,[],2));
end
else
    mi=[]; mod_ind=[];
end
%% compute coherence
% 
% ONLY TMP COMMENTED OUT !!!!!!!
cfg=[];
cfg.output='fourier';
cfg.method='mtmfft';
cfg.taper='dpss';
cfg.foilim=[2 round(data.fsample./2)];%CHANGE
cfg.tapsmofrq=.5;
cfg.keeptrials='yes';
freq=ft_freqanalysis(cfg,data);

cfg=[];
cfg.method='coh';
cfg.complex='complex';
coh=ft_connectivityanalysis(cfg,freq);
cohreal=coh.cohspctrm;
% cfg=[];
% cfg.method='pdc';
% pdc=ft_connectivityanalysis(cfg,freq);
% pdcreal=pdc.pdcspctrm;
cfg=[];
cfg.output='fourier';
cfg.method='mtmfft';
cfg.taper='dpss';
cfg.foilim=[2 25];%CHANGE
cfg.tapsmofrq=.5;
cfg.keeptrials='yes';
freqpsi=ft_freqanalysis(cfg,data);
cfg=[];
cfg.method='psi'; cfg.bandwidth=1;
psi=ft_connectivityanalysis(cfg,freqpsi);
psireal=psi.psispctrm;
nPerms=200;
for perm=1:nPerms
    freqshff=freq;
    for pchan=1:size(freq.fourierspctrm,2)
        freqshff.fourierspctrm(:,pchan,:) = freq.fourierspctrm(randperm(size(freq.fourierspctrm,1)),pchan,:);
    end    
    disp(strcat('Perm',{' '},num2str(perm)))
cfg=[];
cfg.method='coh';
coh=ft_connectivityanalysis(cfg,freqshff);
cohall(:,:,:,perm)=coh.cohspctrm;
% cfg=[];
% cfg.method='pdc';
% pdc=ft_connectivityanalysis(cfg,freqshff);
% pdcall(:,:,:,perm)=pdc.pdcspctrm;
% % 
% cfg=[];
% cfg.method='psi';
% psi=ft_connectivityanalysis(cfg,freqshff);
% psiall(:,:,:,perm)=psi.psispctrm;
% 
end
freq_ax=coh.freq;
% %%
% clear coh psi pdc
cohz=(cohreal-nanmean(cohall,4))./nanstd(cohall,[],4);
% pdc=(pdcreal-nanmean(pdcall,4))./nanstd(pdcall,[],4);
% psi=(psireal-nanmean(psiall,4))./nanstd(psiall,[],4);
% 
% %% compute alpha-gamma power correlation
if data.fsample>=200
[~,alphind(1)]=min(abs(freq.freq-8)); [~,alphind(2)]=min(abs(freq.freq-12));
[~,gammind(1)]=min(abs(freq.freq-70)); [~,gammind(2)]=min(abs(freq.freq-190));%CHANGE
[~,bind(1)]=min(abs(freq.freq-2)); [~,bind(2)]=min(abs(freq.freq-200));%CHANGE

apow=abs(nanmean(freq.fourierspctrm(:,1:(length(data.label)/2),alphind(1):alphind(2)),3));
gpow=abs(nanmean(freq.fourierspctrm(:,1:(length(data.label)/2),gammind(1):gammind(2)),3));
bpow=abs(nanmean(freq.fourierspctrm(:,1:(length(data.label)/2),bind(1):bind(2)),3));


[r, pv]=partialcorr(apow,gpow,bpow);
end
%% hilbert alpha filter find troughs in each channel
for chan = 1:length(lbl)%CHANGE
    disp(strcat('Now processing channel',{' '},num2str(chan),{' '},'out of',{' '},num2str(length(lbl))))
trldf=zeros(length(alph.trial),3);
l_nu=.5;
l_old=2;
clear trldf
pk_counter=1;
for i=1:length(alph.trial)
     [~,pk_inds]=findpeaks(alph.trial{i}(chan,:));
     pk_inds=pk_inds(intersect(find(pk_inds<((l_old-(l_nu/2))*alph.fsample)),find(pk_inds>((l_nu/2)*alph.fsample))));
     if isempty(pk_inds)
        continue
     end
     for ii=1:length(pk_inds)
        
    new_tr_cnt=pk_inds(ii);
    new_tr_cnt=new_tr_cnt+data.sampleinfo(i,1);
    trldf(pk_counter,1)=new_tr_cnt-alph.fsample.*((l_nu/2));
    trldf(pk_counter,2)=new_tr_cnt+alph.fsample.*((l_nu/2));
    trldf(pk_counter,3)=-alph.fsample.*((l_nu/2));%fix
    pk_counter=pk_counter+1;
     end
end
bad_tr_inds=zeros(1,length(trldf));
for trl_ind=1:size(trldf,1)
    if trldf(trl_ind,1)<=0
        bad_tr_inds(trl_ind)=1;
    end
end
trldf(find(bad_tr_inds),:)=[];
  
cfg=[];
cfg.trl=trldf;
alpha_aligned=ft_redefinetrial(cfg,data);

cfg=[]; alpha_avg=ft_timelockanalysis(cfg,alpha_aligned);
alpha_td(chan,:,:)=alpha_avg.avg;
alpha_var(chan,:,:)=alpha_avg.var;
ntrs=length(alpha_avg.cfg.previous.trl);
end