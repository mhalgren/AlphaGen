clear

load('TCData.mat')

%% get chanxtime matrix of how much alpha is greater than power in surrounding bands (as derived from hilbert)

cfg=[]; cfg.bpfreq=[7 13]; cfg.bpfilter='yes'; cfg.bpinstabilityfix='reduce'; alpha=ft_preprocessing(cfg,data);
cfg=[]; cfg.hilbert='abs'; alpha=ft_preprocessing(cfg,alpha);

cfg=[]; cfg.bpfreq=[4 6]; cfg.bpfilter='yes'; cfg.bpinstabilityfix='reduce'; theta=ft_preprocessing(cfg,data);
cfg=[]; cfg.hilbert='abs'; theta=ft_preprocessing(cfg,theta);

cfg=[]; cfg.bpfreq=[14 16]; cfg.bpfilter='yes'; cfg.bpinstabilityfix='reduce'; beta=ft_preprocessing(cfg,data);
cfg=[]; cfg.hilbert='abs'; beta=ft_preprocessing(cfg,beta);

ap=alpha.trial{1};

mp=(beta.trial{1}+theta.trial{1})./2; ma=ap./mp;

mabs=zeros(size(ma)); mabs(ma>=2)=1; ds=diff(mabs,[],2);

mabs2=zeros(size(ma)); mabs2(ma>=1.5)=1; ds2=diff(mabs2,[],2);

bi=1; clear brsti brstp
for k=1:size(ds,1)
    tmpst=ds(k,:); tmpst_i=find(tmpst==1); tmpste=ds2(k,:);
    for kk=tmpst_i
        tmp=tmpste(kk:end); tmpe=find(tmp==-1); 
        if ~isempty(tmpe)
        tmpe=min(tmpe); brsti(bi)=tmpe; brstp(bi)=mean(ma(k,kk:(kk+tmpe))); bi=bi+1;
        end
    end
    disp(k)
end

sum(brsti(brsti>=400))./sum(brsti)