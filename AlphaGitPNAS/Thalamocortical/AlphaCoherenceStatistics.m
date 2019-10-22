%% creates coha
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'Resultsv5.mat'))

    freqb=[7 13];
    [~,lind]=min(abs(freq_ax-freqb(1)));
    [~,hind]=min(abs(freq_ax-freqb(2)));
    
    cohap=squeeze(mean(cohall(:,:,lind:hind,:),3));
    
    cohapmu=squeeze(mean(cohap,3));
    
    cohapstd=squeeze(std(cohap,[],3));
    
    coha=(squeeze(mean(cohreal(:,:,lind:hind),3))-cohapmu)./cohapstd;
    
    save(strcat('L',num2str(k),'coha.mat'),'coha','lbl');
    
    clear
end

%% visualize coha
sind=1;
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'coha.mat'))
    subplot(3,3,sind);
    imagesc(coha), colorbar, set(gca,'XTick',1:length(lbl),'YTick',1:length(lbl),'XTickLabel',lbl,'YTickLabel',lbl)
    sind=sind+1;
end

%% 
sind=1;
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'coha.mat'))
    load(strcat('L',num2str(k),'Resultsv5.mat'),'lbl')
    
    if size(coha,1)>length(lbl)
        
    
    lbl2=lbl;
    for kk=1:length(lbl)
        lbl2{kk}=[lbl{kk} 'HGP'];
    end
    lbl=[lbl lbl2];
    end
    subplot(3,3,sind);
    
    b=(size(coha,1)*(size(coha,1)-1))./2;
    b=.05./b;
    z=-sqrt(2) * erfcinv(b*2);
    cohasig=coha>=-z;
    
    imagesc(cohasig), set(gca,'XTick',1:length(lbl),'YTick',1:length(lbl),'XTickLabel',lbl,'YTickLabel',lbl)
    sind=sind+1;
    clear lbl lbl2
end
suptitle('Coherence Stats')

%% 
clear
sind=1;
ttc=0; ct=0; tc=0; ccs=0; tts=0; cts=0; tcs=0;
figure
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'coha.mat'))
    load(strcat('L',num2str(k),'Resultsv5.mat'),'lbl')
    load(strcat('L',num2str(k),'ThChs.mat'))
    
    if size(coha,1)>length(lbl)
    if k==1
        coha=coha(2:end,2:end); ThChs=ThChs-1; lbl=lbl(2:end);
    end
    CtxChs=1:length(lbl); CtxChs(ThChs)=[];
    
    lbl=lbl(:);
    
    lbl = [lbl(ThChs); lbl(CtxChs)];
        
    lbl2=lbl;
    for kk=1:length(lbl)
        lbl2{kk}=[lbl{kk} 'HGP'];
    end
    lbl=[lbl lbl2]; lbl=lbl(:);
    
    
    b=length(ThChs)*length(lbl2);%b=(size(coha,1)*(size(coha,1)-1))./2;
    b=.05./b;
    z=-sqrt(2) * erfcinv(b*2);
    cohasig=coha>=-z;
    
    tmp=cohasig(ThChs,ThChs+length(lbl2));
    %%
    tctmp=cohasig(ThChs,CtxChs+length(lbl2)); cttmp=cohasig(CtxChs,ThChs+length(lbl2)); tttmp=cohasig(ThChs,ThChs+length(lbl2)); cctmp=cohasig(CtxChs,CtxChs+length(lbl2));
    tc=tc+sum(tctmp(:));
    ct=ct+sum(cttmp(:));
    
    tttmp=sum(tttmp,2); tts=tts+length(find(tttmp>0)); clear tttmp
    cctmp=sum(cctmp,2); ccs=ccs+length(find(cctmp>0)); clear cctmp
    tctmp=sum(tctmp,2); tcs=tcs+length(find(tctmp>0)); clear tctmp
    cttmp=sum(cttmp,2); cts=cts+length(find(cttmp>0));
    %%
    ttc=ttc+length(ThChs).^2;
    
    coha=coha(:,(length(lbl2)+1):end);
    coha=coha(:,[ThChs CtxChs]);
    
    cohasig=cohasig(:,(length(lbl2)+1):end);
    cohasig=cohasig(:,[ThChs CtxChs]);
    
    subplot(2,3,sind), sind=sind+1;
    imagesc(coha(ThChs,:),[0 4]), set(gca,'YTick',1:length(ThChs),'XTick',[],'XTickLabel',[],'YTickLabel',lbl(1:length(ThChs))), hold on
    %imcontour(gca,1:size(cohasig,2),1:length(ThChs),cohasig(ThChs,:),1,'linecolor','k','linew',2), box off, axis equal
    imcontour(cohasig(ThChs,:),1), box off
    title(num2str(sum(cttmp(:))))
    %title(num2str(sum(tmp(:))))
    end
    clear lbl lbl2
end



clear
figure
sind=1; tc=0; ct=0; cc=0; ccc=0; ccs=0; tts=0; cts=0; tcs=0;
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'Resultsv5.mat'),'mod_ind','lbl')
    load(strcat('L',num2str(k),'ThChs.mat'))
    
    if k==1
        mod_ind=mod_ind(2:end,2:end); ThChs=ThChs-1; lbl=lbl(2:end);
    end
    
    CtxChs=1:length(lbl); CtxChs(ThChs)=[];
    
    lbl=lbl(:);
    
    lbl = [lbl(ThChs); lbl(CtxChs)];        
    if ~isempty(mod_ind)
    
    subplot(2,3,sind);
    sind=sind+1;
    b=(length(ThChs)*(size(mod_ind,1)));%b=(size(mod_ind,1)*(size(mod_ind,1)));
    b=.05./b;
    z=-sqrt(2) * erfcinv(b*2);
    cohasig=mod_ind>=-z;
    
    mod_ind=mod_ind(:,[ThChs CtxChs]);
    
    imagesc(mod_ind(ThChs,:),[0 4]), set(gca,'XTick',[],'YTick',1:length(ThChs),'XTickLabel',[],'YTickLabel',lbl(1:length(ThChs))), box off, hold on
    imcontour(cohasig(ThChs,:),1), box off
    tmp=cohasig(ThChs,ThChs);
    %%
    tctmp=cohasig(ThChs,CtxChs); cttmp=cohasig(CtxChs,ThChs); cctmp=cohasig(CtxChs,CtxChs); tttmp=cohasig(ThChs,ThChs);
    tc=tc+sum(tctmp(:));
    ct=ct+sum(cttmp(:));
    cc=cc+sum(cctmp(:));
    
    ccc=ccc+length(CtxChs);
    
    tttmp=sum(tttmp,2); tts=tts+length(find(tttmp>0)); disp(tts), clear tttmp
    cctmp=sum(cctmp,2); ccs=ccs+length(find(cctmp>0)); clear cctmp
    tctmp=sum(tctmp,2); tcs=tcs+length(find(tctmp>0)); clear tctmp
    cttmp=sum(cttmp,2); cts=cts+length(find(cttmp>0));    
    
    %%
    title(num2str(sum(cttmp(:))))
    %title(num2str(sum(tmp(:))))
    end
end


sind=1;
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'Resultsv5.mat'),'mod_ind','lbl')
    subplot(3,3,sind);
    
    imagesc(mod_ind,[0 5]), set(gca,'XTick',1:length(lbl),'YTick',1:length(lbl),'XTickLabel',lbl,'YTickLabel',lbl)
    sind=sind+1;
end