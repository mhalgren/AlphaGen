function [SRmu,SRstd] = SharpnessRatio2(data)

d=fieldtrip2mat(data);

fs=data.fsample;

cfg=[]; cfg.bpfilter=[7 13]; cfg.bpinstabilityfix='reduce'; data=ft_preprocessing(cfg,data);

%cfg=[]; cfg.hilbert='abs'; pow=ft_preprocessing(cfg,data);

for ch=1:size(data.trial{1},1) %loop through channels
    sh=[];
%     p=pow.trial{1}(ch,:);
%     [~,p]=sort(p);
%     p=p(round(length(p)*.9):end);
    
    dat=d(ch,:);
    zc=crossing(dat);
    %zc=intersect(zc,p);
    ts=round(.005/(1/fs)); %find number of sample to go out from extrema to calculate sharpness

for k=1:(length(zc)-1) %loop thru pos/negative segments
    x=dat((zc(k)+1):zc(k+1));
    if abs(max(x))>abs(min(x))
        [~,ex]=max(x);
    else
        [~,ex]=min(x);
    end
    ex=ex+zc(k)-1;
    if (ex-ts)>=1 && (ex+ts)<=length(dat)
        if isempty(sh)
            sh=abs((dat(ex)-dat(ex-ts)) + (dat(ex)-dat(ex+ts)))/2;
        else
            sh=[sh abs((dat(ex)-dat(ex-ts)) + (dat(ex)-dat(ex+ts)))/2];
        end;
    end
end

pk=sh(1:2:end); tr=sh(2:2:end); %arbitrary division into peak/trough
len = min([length(pk) length(tr)]);
pk = pk(1:len); tr = tr(1:len);

[~,SRtmp] = max([nanmean(pk)./nanmean(tr) nanmean(tr)./nanmean(pk)]);
if SRtmp==1
SRmu(ch) =  mean(pk./tr);
SRstd(ch) = std(pk./tr);
elseif SRtmp==2
SRmu(ch) =  mean(tr./pk);
SRstd(ch) = std(tr./pk);
end
disp(ch)
end
end