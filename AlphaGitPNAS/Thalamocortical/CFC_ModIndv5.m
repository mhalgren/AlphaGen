function[mi,mod_ind] = CFC_ModIndv5(data,frqs)

nChans=length(data.label);

mi=zeros(length(frqs),nChans,nChans,18,length(frqs));

mod_ind=zeros(length(frqs),nChans,nChans,length(frqs));

cfg=[];
cfg.bsfilter='yes';
cfg.bsfreq=[48 52;98 102; 148 152;198 202;248 252;398 402;448 452;498 502];
cfg.bsinstabilityfix='reduce';
data=ft_preprocessing(cfg,data);

allamplitude=zeros(nChans,length(data.trial{1}(:,1:10:end)),length(frqs));

for kk=1:length(frqs)
    disp(kk)
cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=frqs(kk,:);
cfg.bpinstabilityfix='reduce';
amp_data=ft_preprocessing(cfg,data);
cfg=[];
cfg.hilbert='abs';
amp_data=ft_preprocessing(cfg,amp_data);
tmp=fieldtrip2mat(amp_data);
allamplitude(:,:,kk)=squeeze(tmp(:,1:10:end));
end

shuffpoint=zeros(1,100);
for z=1:100
    shuffpoint(z)=randi(round(length(allamplitude).*.6))+round(length(allamplitude).*.2);
end

for k=1:length(frqs)
%% only use alpha trials
cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=frqs(k,:);
cfg.bpinstabilityfix='reduce';
alph=ft_preprocessing(cfg,data);
cfg=[];
cfg.hilbert='angle';
phs_data=ft_preprocessing(cfg,alph);
allphase=fieldtrip2mat(phs_data);
allphase=allphase(:,1:10:end);
%%
%filter in phase-freq
for kk=1:length(frqs)
for chan = 1:length(data.label)
disp(strcat('Now processing channel',{' '},num2str(chan),{' '},'out of',{' '},num2str(length(data.label))))
%calculate MI
    phase=squeeze(allphase(chan,:));

    real_pac=average_envelope_versus_phase_v4(allamplitude,phase);
    tmp=zeros(size(allamplitude,1),18,size(allamplitude,3),100);
    
    for nperm=1:100
        disp(nperm)
        shuffphase=squeeze(phase([shuffpoint(nperm):end 1:(shuffpoint(nperm)-1)]));
        tmp(:,:,:,nperm)=average_envelope_versus_phase_v4(allamplitude,shuffphase);
        %normalize tmp
    end

    mi(k,chan,:,:,:)=squeeze((real_pac-mean(tmp,4))./std(tmp,[],4));
    tmp=tmp./repmat(sum(tmp,4),[1 1 1 size(tmp,4)]);
    real_pac=real_pac./repmat(sum(real_pac,4),[1 1 1 size(real_pac,4)]);
    real_pac2=squeeze((log(18)+sum(real_pac.*log(real_pac),2))./log(18));
    tmp2=squeeze((log(18)+sum(tmp.*log(tmp),2))./log(18));
    mod_ind(k,chan,:,:)=squeeze((real_pac2-mean(tmp2,3))./std(tmp2,[],3));
end
disp(num2str(kk))
end
disp(num2str(k))
end