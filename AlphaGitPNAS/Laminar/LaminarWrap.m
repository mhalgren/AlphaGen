clear
wdir='';
addpath(genpath(strcat(wdir,'Laminar Oscillations/Miscellaneous Subfunctions/fieldtrip-20150121')))
cd(strcat(wdir,'AlphaLaminar'))
addpath(genpath('Scripts'))
addpath(genpath('PGData'))
addpath(genpath('MUA_PreProcessed'))

%filename
fnms={'NYEOEC1#1@1e.cnt','02a_FF.cnt','HungLam.mat',};

snms={};

epoched=[0 0 0];

%ten_sec_good - 1 if BAD_TRIALS ARE ACTUALLY GOOD TRIALS
ten_sec_good=[0 0 0 1 1 0 0 2];

is_hungarian=[0 0 0 0 0 1 0 0];

hasmua=[1 0 1];% ADD MUA FOR NY09!!!!!!
%output...
%whole pg
%whole csd
%trialdef
%alpha on peak
%evoked tfr

csd_method=0; 

for sub_ind=1:3
    disp(strcat('Subject',{' '},num2str(sub_ind)))
    snm=snms{sub_ind};
    load(strcat('ChLbls/',snm,'ChLbls.mat'))
    if sub_ind==6
        chlbl=chlbl(1:17);
    end
    if hasmua(sub_ind)
    load(strcat(snm,'_RawMUA.mat'))
    %mua.trial{1}=zscore(mua.trial{1},[],2);
    %h=fspecial('gaussian',[1 50],20);
    %mua.trial{1}=imfilter(mua.trial{1},h);
    mua.trial{1}=abs(mua.trial{1});
    end
    %segment epochs
    cfg=[];
    if ~is_hungarian(sub_ind) && sub_ind~=8
    cfg.dataset=fnms{sub_ind};
    elseif sub_ind==8 %&& sub_ind==7
            load(fnms{sub_ind})
    else
        load('HungLam.mat')
    end
    if epoched(sub_ind)~=1
        cfg.continuous='yes';
    else
        cfg.continuous='no';
    end
    cfg.bsfilter='yes';
    if is_hungarian(sub_ind)
        cfg.bsfreq=[48 52; 98 102; 148 152; 198 202; 248 252];
    else
    cfg.bsfreq=[58 62; 118 122; 178 182; 238 242];
    end
    cfg.bsinstabilityfix='reduce';
    cfg.hpfilter='yes';
    if sub_ind~=5
    cfg.hpfreq=.5;
    else
        cfg.hpfreq=3.5;
    end
    cfg.hpinstabilityfix='reduce';
    cfg.channel=chlbl;
    if sub_ind==5
        cfg.channel={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19'};
    end
    if ~is_hungarian(sub_ind) && sub_ind~=8 %&& sub_ind~=7
    pgdata=ft_preprocessing(cfg);
    else
        pgdata=ft_preprocessing(cfg,data);
    end
    if sub_ind==7
        pgdata.trial{1}(1,:)=pgdata.trial{1}(3,:).*(1./3);
        pgdata.trial{1}(2,:)=pgdata.trial{1}(3,:).*(2./3);

        pgdata.trial{1}(7,:)=((2./3).*pgdata.trial{1}(6,:)) + ((1./3).*pgdata.trial{1}(9,:));
        pgdata.trial{1}(8,:)=((1./3).*pgdata.trial{1}(6,:)) + ((2./3).*pgdata.trial{1}(9,:));

        pgdata.trial{1}(11,:)=((4./5).*pgdata.trial{1}(10,:)) + ((1./5).*pgdata.trial{1}(15,:));
        pgdata.trial{1}(12,:)=((3./5).*pgdata.trial{1}(10,:)) + ((2./5).*pgdata.trial{1}(15,:));
        pgdata.trial{1}(13,:)=((2./5).*pgdata.trial{1}(10,:)) + ((3./5).*pgdata.trial{1}(15,:));
        pgdata.trial{1}(14,:)=((1./5).*pgdata.trial{1}(10,:)) + ((4./5).*pgdata.trial{1}(15,:));

        pgdata.trial{1}(20,:)=((2./3).*pgdata.trial{1}(19,:)) + ((1./3).*pgdata.trial{1}(22,:));
        pgdata.trial{1}(21,:)=((1./3).*pgdata.trial{1}(19,:)) + ((2./3).*pgdata.trial{1}(22,:));
    end
    raw_upperlim=pgdata.sampleinfo(end,2);
    if sub_ind==6
        load(strcat(wdir,'AlphaLaminar/Grid/ecogdata.mat'))
        ecog_dat.sampleinfo=pgdata.sampleinfo;
        cfg=[];
        cfg.continuous='yes';
        cfg.bsfreq=[48 52; 98 102; 148 152; 198 202; 248 252];
        cfg.bsinstabilityfix='reduce';
        cfg.hpfilter='yes';
        cfg.hpfreq=.5;
        cfg.hpinstabilityfix='reduce';
        ecog_dat=ft_preprocessing(cfg,ecog_dat);
    end
    switch epoched(sub_ind)
        case 0
            cfg=[]; cfg.length=5; cfg.overlap=0;
            pgdata=ft_redefinetrial(cfg,pgdata);
            if sub_ind==6
            ecog_dat=ft_redefinetrial(cfg,ecog_dat);
            end
            if hasmua(sub_ind)
            mua=ft_redefinetrial(cfg,mua);
            end
        case 2
            cfg=[]; cfg.length=10; cfg.overlap=0;
            pgdata=ft_redefinetrial(cfg,pgdata);
            if hasmua(sub_ind)
            mua=ft_redefinetrial(cfg,mua);
            end
        case 1
            %
    end
    load(strcat('BadTrials/',snm,'BadTrials.mat'))
    
    if ten_sec_good(sub_ind)==0
        
    %reject bad trials
    cfg=[]; cfg.trials=1:length(pgdata.trial);
    if sub_ind==7
        bad_trials=bad_trials(1:11);
    end
    cfg.trials(bad_trials)=[];
    pgdata=ft_redefinetrial(cfg,pgdata);
    if sub_ind==6
    ecog_dat=ft_redefinetrial(cfg,ecog_dat);
    end
    if hasmua(sub_ind)
    mua=ft_redefinetrial(cfg,mua);
    muaallmu=fieldtrip2mat(mua);
    muaallmu=mean(muaallmu,2);
    end
    elseif ten_sec_good(sub_ind)==1
    cfg=[]; cfg.trials=bad_trials;
    pgdata=ft_redefinetrial(cfg,pgdata);
    if hasmua(sub_ind)
    mua=ft_redefinetrial(cfg,mua);
    muaallmu=fieldtrip2mat(mua);
    muaallmu=mean(muaallmu,2);
    end
    else
        if hasmua(sub_ind)
            muaallmu=fieldtrip2mat(mua);
            muaallmu=mean(muaallmu,2);
        end
    end 

    if sub_ind==3
        for k=1:length(pgdata.trial)

        pgdata.trial{k}(14,:)=(1/2).*(pgdata.trial{k}(13,:)+pgdata.trial{k}(15,:));
        end
    end
    if sub_ind==5
        for k=1:length(pgdata.trial)
        pgdata.trial{k}(1,:)=.5.*pgdata.trial{k}(2,:);
        pgdata.trial{k}(3,:)=.5.*(pgdata.trial{k}(2,:)+pgdata.trial{1}(4,:));
        
        pgdata.trial{k}(12,:)=.5.*(pgdata.trial{k}(11,:)+pgdata.trial{1}(13,:));
        
        end
    end
    %convert to csd and ecd
    csdata=pgdata;
    ecdata=pgdata;
    csdata.label=pgdata.label(2:end-1);
    ecdata.label={'ecd'};
    csdata.hdr.nChans=length(csdata.label);
    ecdata.hdr.nChans=1;
    csdata.hdr.label=csdata.label;
    ecdata.hdr.label=ecdata.label;
    if sub_ind~=5
        hd=.15;
    else
        hd=.175;
    end
    for k=1:length(csdata.trial)
        if csd_method
            [ecdata.trial{k},csdata.trial{k},~]=ts_calc_laminar_ecd(pgdata.trial{k},pgdata.fsample);%sergei's method
        else
            csdata.trial{k}=pg2csdv3(pgdata.trial{k},hd);%istvan's method
            [ecdata.trial{k},~,~]=ts_calc_laminar_ecd(pgdata.trial{k},pgdata.fsample);
        end
    end
    %% only use trials w/ alpha bump
cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[7 13];
cfg.bpinstabilityfix='reduce';
alph=ft_preprocessing(cfg,csdata);

cfg=[];
cfg.hilbert='abs';
alph_amp=ft_preprocessing(cfg,alph);

a=fieldtrip2mat(alph_amp);

cfg=[];
%cfg.output='power';
cfg.taper='dpss';
cfg.keeptrials='yes';
cfg.method='mtmfft';
cfg.tapsmofrq=.2;
cfg.foilim=[1 45];
frqalph=ft_freqanalysis(cfg,csdata);

%% select 2s trials w/ 20% greatest alpha
cfg=[];
cfg.length=2;
cfg.overlap=0;
csdata=ft_redefinetrial(cfg,csdata);
ecdata=ft_redefinetrial(cfg,ecdata);
pgdata=ft_redefinetrial(cfg,pgdata);
if sub_ind==6
ecog_dat=ft_redefinetrial(cfg,ecog_dat);
end
alph=ft_redefinetrial(cfg,alph);
if hasmua(sub_ind)
    mua=ft_redefinetrial(cfg,mua);
end

cfg=[];
cfg.hilbert='abs';
alph_amp=ft_preprocessing(cfg,alph);

cfg=[];
cfg.hilbert='angle';
alph_phs=ft_preprocessing(cfg,alph);

a=squeeze(mean(squeeze(fieldtrip2mat_epochs(alph_amp)),2)); a=mean(a);
[~,a_chan]=max(mean(squeeze(mean(squeeze(fieldtrip2mat_epochs(alph_amp)),2)),2));

if sub_ind==3 
a_chan=3;
end

if ~isequal(snm,'')
[~,ind]=sort(a); ind=ind(round(length(ind).*.8):end);
    % reject non alpha epochs
    cfg=[];
    cfg.trials=ind;
    ecdata=ft_redefinetrial(cfg,ecdata);
    csdata=ft_redefinetrial(cfg,csdata);
    pgdata=ft_redefinetrial(cfg,pgdata);
    if sub_ind==6
    ecog_dat=ft_redefinetrial(cfg,ecog_dat);
    end
    alph=ft_redefinetrial(cfg,alph);
    alph_amp=ft_redefinetrial(cfg,alph_amp);
    alph_phs=ft_redefinetrial(cfg,alph_phs);
    if hasmua(sub_ind)
    mua=ft_redefinetrial(cfg,mua);
    muamu=fieldtrip2mat(mua);
    muamu=mean(muamu,2);
    end
end

    %%
    %hilbert alpha filter find troughs in each channel
    cfg=[];
    cfg.bpfilter='yes';
    cfg.bpfreq=[7 13];
    cfg.bpinstabilityfix='reduce';
    ecdalph=ft_preprocessing(cfg,ecdata);
    cfg=[];
    cfg.hilbert='angle';
    ecdalph=ft_preprocessing(cfg,ecdalph);
%     trldf=zeros(length(alph.trial),3);
    clear trldf
    pk_counter=1;
    for i=1:length(alph.trial)
        [~,pk_inds]=findpeaks(alph.trial{i}(a_chan,:));
        %[~,pk_inds]=findpeaks(ecdalph.trial{i});
        pk_inds=pk_inds(intersect(find(pk_inds<(1.75*ecdata.fsample)),find(pk_inds>(.25*ecdata.fsample))));
        if isempty(pk_inds)
            continue
        end
        for ii=1:length(pk_inds)
        new_tr_cnt=pk_inds(ii);
        new_tr_cnt=new_tr_cnt+alph.sampleinfo(i,1);
        trldf(pk_counter,1)=new_tr_cnt-alph.fsample.*(.25);
        trldf(pk_counter,2)=new_tr_cnt+alph.fsample.*(.25);
        trldf(pk_counter,3)=-alph.fsample.*(.25);%fix
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
    
    if isempty(trldf)
    disp(strcat('No alpha in subject',{' '},snm,', continuing'))
    continue
    end
    
    cfg=[];
    cfg.trl=trldf;
    alpha_aligned=ft_redefinetrial(cfg,csdata);
    ecdata_aligned=ft_redefinetrial(cfg,ecdata);
    pgdata_aligned=ft_redefinetrial(cfg,pgdata);
    if sub_ind==6
    ecog_aligned=ft_redefinetrial(cfg,ecog_dat);
    end
    cfg=[]; alpha_avg=ft_timelockanalysis(cfg,alpha_aligned);
    cfg=[]; ecdata_avg=ft_timelockanalysis(cfg,ecdata_aligned);
    if sub_ind==6
    cfg=[]; ecog_avg=ft_timelockanalysis(cfg,ecog_aligned);
    end
    %% find fourier spectrm of trough-locked epochs
    cfg=[];
    cfg.output='fourier';
    cfg.method='mtmfft';
    cfg.taper='dpss';
    cfg.foilim=[3 200];
    cfg.tapsmofrq=2;
    cfg.keeptrials='yes';
    freqecd=ft_freqanalysis(cfg,ecdata_aligned);
    freqalph=ft_freqanalysis(cfg,alpha_aligned);
    %% find evoked power

    cfg=[];
    cfg.method='wavelet';
    cfg.foilim=[15 300];
    cfg.toi=[-.25:.01:.25];
    alpha_tfr=ft_freqanalysis(cfg,pgdata_aligned);
    
    %% find coherence
    
    cfg=[];
    cfg.bpfilter='yes';
    cfg.bpfreq=[70 190];
    hgp=ft_preprocessing(cfg,pgdata);
    cfg=[];
    cfg.hilbert='abs';
    hgp=ft_preprocessing(cfg,hgp);
    
    hmu=fieldtrip2mat(hgp);
    hmu=nanmean(hmu,2);
    
    for k=1:length(hgp.label)
        hgp.label{k}=strcat('h',hgp.label{k});
        hgp.cfg.channel{k}=hgp.label{k};
    end
    cfg=[];
    concat_data=ft_appenddata(cfg,csdata,hgp);
    
    cfg=[];
    cfg.trl=trldf;
    hgp_avg=ft_redefinetrial(cfg,hgp);
    cfg=[];
    hgp_avg=ft_timelockanalysis(cfg,hgp_avg);
    
    if hasmua(sub_ind)
    cfg=[];
    cfg.trl=trldf;
    mua_aligned=ft_redefinetrial(cfg,mua);
    mua_avg=ft_timelockanalysis(cfg,mua_aligned);
    save(strcat('LoopOut10/',snm,'_AvgMUA.mat'),'mua_avg')

    for k=1:length(mua.label)
        mua.label{k}=strcat('m',mua.label{k});
        mua.cfg.channel{k}=mua.label{k};
    end
    
    cfg=[];
    concat_data=ft_appenddata(cfg,concat_data,mua);
%     cfg=[];
%     cfg.length=1;
%     cfg.overlap=.5;
%     concat_data=ft_redefinetrial(cfg,concat_data);
%else    
%     cfg=[];
%     cfg.length=1;
%     cfg.overlap=.5;
%     concat_data=ft_redefinetrial(cfg,concat_data);
    end
       
    cfg=[];
    cfg.output='fourier';
    cfg.method='mtmfft';
    cfg.taper='dpss';
    cfg.foilim=[2 200];
    cfg.tapsmofrq=.5;
    cfg.keeptrials='yes';
    freq=ft_freqanalysis(cfg,concat_data);
    freqecdlng=ft_freqanalysis(cfg,ecdata);
    cfg=[];
    cfg.method='coh';
    coh=ft_connectivityanalysis(cfg,freq);
    cohreal=coh.cohspctrm;
    nPerms=200;%change to 200???
    for perm=1:nPerms
    freqshff=freq;
    for pchan=1:size(freq.fourierspctrm,2)
    freqshff.fourierspctrm(:,pchan,:) = freq.fourierspctrm(randperm(size(freq.fourierspctrm,1)),pchan,:);
    end
    
    disp(strcat('Perm',{' '},num2str(perm)))
    cfg=[];
    cfg.method='coh';
    cohp=ft_connectivityanalysis(cfg,freqshff);
    cohall(:,:,:,perm)=cohp.cohspctrm;
    end
    cohz=(cohreal-nanmean(cohall,4))./nanstd(cohall,[],4);
    
    PAC=NaN(length(hgp.label),length(csdata.label),37);
    h=fieldtrip2mat(hgp); pa=fieldtrip2mat(alph_phs); pae=fieldtrip2mat(ecdalph);
                for ch_ind=1:(length(csdata.label)+1)
                    if ch_ind<=length(csdata.label)
                        p=pa(ch_ind,:);
                    else
                        p=pae;
                    end
              disp(ch_ind)
      p = p * 360 / (2.0*pi);              
      j=1;
      for k=-180:10:180
          inds = intersect(find(p >= k), find(p < k+10));  
          PAC(:,j) = squeeze(mean(h(:,inds),2));  
          j=j+1;
      end
      nPerms=200;
      PACz=NaN(length(hgp.label),37,nPerms);
      for perm=1:nPerms
          shuffpoint=randi(round(length(p).*.6))+round(length(p).*.2);
          p_shuff=p([shuffpoint:end 1:(shuffpoint-1)]);
          j=1;
      for k=-180:10:180
          inds = find(p_shuff >= k & p_shuff < k+10);  %Find phase values in a 10 degree band.
          PACz(:,j,perm) = squeeze(mean(h(:,inds),2));       %Compute mean amplitude at these phases.
          j=j+1;
      end
      %disp(perm)
      end
      
      PAC=PAC(:,1:36);
      PACz=squeeze(PACz(:,1:36,:));
      
      %% calculate Modulation Index
      
      PAC=PAC./repmat(sum(PAC,2),[1 36]);
      PACz=PACz./repmat(sum(PACz,2),[1 36]);
      
      ModInd=squeeze(sum(PAC.*log(PAC),2));
      ModIndNull=squeeze(sum(PACz.*log(PACz),2));
      
      ModIndZ=(ModInd-mean(ModIndNull,2))./std(ModIndNull,[],2);
      
      PACa(ch_ind,:,:)=PAC;
      PACza(ch_ind,:,:,:)=PACz;
      ModInda(ch_ind,:)=ModInd;
      ModIndZa(ch_ind,:)=ModIndZ;
    end
      %% Do MUA Mod Ind
      if hasmua(sub_ind)
                    PACm=NaN(length(mua.label),37);
h=fieldtrip2mat(mua);
                for ch_ind=1:(length(csdata.label)+1)
                    if ch_ind<=length(csdata.label)
                        p=pa(ch_ind,:);
                    else
                        p=pae;
                    end
                    disp(ch_ind)

      p = p * 360.0 / (2.0*pi);                 %Use phase in degrees.
      j=1;
      for k=-180:10:180
          inds = find(p >= k & p < k+10);  %Find phase values in a 10 degree band.
          PACm(:,j) = squeeze(mean(h(:,inds),2));       %Compute mean amplitude at these phases.
          j=j+1;
      end
      nPerms=100;
      PACzm=NaN(length(mua.label),37,nPerms);
      for perm=1:nPerms
          shuffpoint=randi(round(length(p).*.6))+round(length(p).*.2);
          p_shuff=p([shuffpoint:end 1:(shuffpoint-1)]);
          j=1;
      for k=-180:10:180
          inds = find(p_shuff >= k & p_shuff < k+10);  %Find phase values in a 10 degree band.
          PACzm(:,j,perm) = squeeze(mean(h(:,inds),2));       %Compute mean amplitude at these phases.
          j=j+1;
      end
      end
      
      PACm=PACm(:,1:36);
      PACzm=squeeze(PACzm(:,1:36,:));
      
      %% calculate Modulation Index
      
      PACm=PACm./repmat(sum(PACm,2),[1 36]);
      PACzm=PACzm./repmat(sum(PACzm,2),[1 36]);
      
      ModIndm=squeeze(sum(PACm.*log(PACm),2));
      ModIndNullm=squeeze(sum(PACzm.*log(PACzm),2));
      
      ModIndZm=(ModIndm-mean(ModIndNullm,2))./std(ModIndNullm,[],2);
      
      PACma(ch_ind,:,:)=PACm;
      PACzma(ch_ind,:,:,:)=PACzm;
      ModIndma(ch_ind,:)=ModIndm;
      ModIndZma(ch_ind,:)=ModIndZm;  
      end
    end
    
    %% save everything
    %save(strcat('LoopOut4/',snm,'_PGData.mat'),'pgdata')
    %save(strcat('LoopOut4/',snm,'_CSData.mat'),'csdata')
    save(strcat('LoopOut10/',snm,'_ECDAvg.mat'),'ecdata_avg')    
    save(strcat('LoopOut10/',snm,'_TrlDf.mat'),'trldf')
%   save(strcat('LoopOut/',snm,'_AlignedAlpha.mat'),'alpha_aligned')
    save(strcat('LoopOut10/',snm,'_AvgAlpha.mat'),'alpha_avg','hgp_avg')
    save(strcat('LoopOut10/',snm,'_AlphaCoh.mat'),'coh','cohz','cohall','-v7.3')
    save(strcat('LoopOut10/',snm,'_AlphaFRQ.mat'),'freq')
    save(strcat('LoopOut10/',snm,'_AlphaVertTW.mat'),'ph')
    save(strcat('LoopOut10/',snm,'_AlphaFreqNu.mat'),'frqalph')
    save(strcat('LoopOut10/',snm,'_AlphaModInd.mat'),'PACa','PACza','ModInda','ModIndZa','hmu')
    if hasmua(sub_ind)
        save(strcat('LoopOut10/',snm,'_AlphaModIndMUA.mat'),'PACma','PACzma','ModIndma','ModIndZma')
    end
    save(strcat('LoopOut10/',snm,'_AlphaFreqMisc.mat'),'freqecd','freqalph','freqecdlng')
    if hasmua(sub_ind)
        save(strcat('LoopOut10/',snm,'_AlphaMUABl.mat'),'muamu','muaallmu')
    end
    save(strcat('LoopOut10/',snm,'_AlphaTFR.mat'),'alpha_tfr')
    if sub_ind==6
        save('LoopOut10/Pt2EcogStuff.mat','ecog_dat','ecog_aligned','ecog_avg')
    disp(strcat('DONE W/ SUBJECT',{' '},num2str(sub_ind)))
    
    end
    clear PACma PACzma ModIndma ModIndZma PACa PACza ModInda ModIndZa cohall
end