function frq = evokedpwr(data,trl,tmp)
if ~isempty(trl)
cfg=[];
cfg.trl=trl;
data=ft_redefinetrial(cfg,data);
end

cfg=[]; cfg.length=3; data=ft_redefinetrial(cfg,data);

if ~isempty(tmp)
art_smpl=tmp.artfctdef.visual.artifact;
bad_trials=zeros(1,length(art_smpl));
for kk=1:size(art_smpl,1);
    smp=round((art_smpl(kk,1)+art_smpl(kk,2))./2);
    for tind=1:length(data.sampleinfo)
        if smp<=data.sampleinfo(tind,2) && smp>=data.sampleinfo(tind,1)
        bad_trials(kk)=tind;
        break
        end
    end
end

bad_trials(find(bad_trials==0))=[];

else
    bad_trials=[];
end
cfg=[];
trs=1:length(data.trial);
trs(bad_trials)=[];
cfg.trials=trs;
data=ft_redefinetrial(cfg,data);

cfg=[]; cfg.method='mtmfft'; cfg.tapsmofrq=1./3; cfg.foilim=[2 45];frq=ft_freqanalysis(cfg,data);
end