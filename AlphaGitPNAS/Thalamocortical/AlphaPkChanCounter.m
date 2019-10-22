clear
wdir='';
addpath(genpath([wdir 'fieldtrip-20150121/']))
cd(strcat(wdir,'AlphaLaminar/TC'))
addpath(genpath(strcat(wdir,'AlphaLaminar/Scripts/')))
addpath(strcat(wdir,'AlphaLaminar/TC/TCPrelimData'))
lbltot={}; ctxtot=0; thtot=0; ctxa=0; tha=0; allch=[]; pkfrq=nan(9,30); si=1;
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'Resultsv5.mat'),'lbl','freq','freq_ax')
    load(strcat('L',num2str(k),'ThChs.mat'))
    
    CtxChs=1:length(lbl); CtxChs(ThChs)=[];
    
    ctxtot=ctxtot+length(CtxChs); thtot=thtot+length(ThChs);
    
    x=freq.fourierspctrm; x=x.*conj(x); x=squeeze(nanmean(x));
    
    gch=[]; pi=1;
    for kk=1:length(lbl)
        f=x(kk,:); locsi=peakfinder(f,((max(f)-min(f))/5)); locs=freq_ax(locsi);
        if ~isempty(find(abs(locs-10)<=3))
            gch = [gch kk];
            gi=find(abs(locs-10)<=3); locs=locs(gi); locsi=locsi(gi);
            [~,maxpk]=max(f(locsi));
            pkfrq(si,pi)=locs(maxpk); pi=pi+1;
        end
    end
    
    ch=zeros(1,length(lbl)); ch(gch)=1;
    
    allch=[allch ch];
    
    lbltot=[lbltot; lbl(:)];
    
    ctxa=ctxa+length(intersect(CtxChs,gch)); tha=tha+length(intersect(ThChs,gch));
    
    si=si+1;
    
    disp(strcat('Done W/',{' '},num2str(k)))
end