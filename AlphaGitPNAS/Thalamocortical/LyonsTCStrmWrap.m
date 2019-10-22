clear
wdir='';
addpath(genpath([wdir 'fieldtrip-20150121/']))
cd(strcat(wdir,'AlphaLaminar/TC'))
addpath(genpath(strcat(wdir,'AlphaLaminar/Scripts/')))
addpath(strcat(wdir,'AlphaLaminar/TC/TCPrelimData'))
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'.mat'))
    if k<=15
        load(strcat('L',num2str(k),'Chs.mat'))
    else
        lbl=data.label;
    end
    load(strcat('L',num2str(k),'BadTrials.mat'))
    load(strcat('L',num2str(k),'ThChs.mat'))
    art_smpl=tmp.artfctdef.visual.artifact;
    bad_trials=zeros(1,length(art_smpl));
    for kk=1:size(art_smpl,1);
        smp=round((art_smpl(kk,1)+art_smpl(kk,2))./2);
        bad_trials(kk)=ceil(smp./(4*data.fsample));
    end
    
    cfg=[]; cfg.channel=lbl; data=ft_preprocessing(cfg,data);
    [alpha_td,lbl,freq,coh,cohz,cohreal,cohall,pdc,psi,freq_ax,r,pv,alpha_var,ntrs,mi,mod_ind,SR] = TC_Strmv4(data,bad_trials,ThChs);    
    save(strcat('L',num2str(k),'Resultsv11.mat'),'alpha_td','lbl','freq','coh','cohz','cohreal','cohall','pdc','psi','freq_ax','r','pv','alpha_var','ntrs','mi','mod_ind','SR','-v7.3')
    disp(strcat('Done W/',{' '},num2str(k)))
    clear
end
