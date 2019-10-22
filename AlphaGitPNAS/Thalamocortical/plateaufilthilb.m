function filt = plateaufilthilb(data,Fs,Fbp,tw)    
filt=data;
for k=1:length(data.trial)
    dat=data.trial{k};
ax = linspace(0, Fs, size(dat,2)); % frequency coefficients
    fl = nearest(ax, min(Fbp))-1; % low cut-off frequency
    fll=round((1-tw)*fl); 
    fh = nearest(ax, max(Fbp))+1; % high cut-off frequency
    fhh=round((1+tw)*fh);
    fc=ones(1,length(ax));
    fc(1:fll)   = 0; % perform low cut-off
    fc(fll:fl) = linspace(0,1,length(fll:fl));
    fc(fhh:end) = 0; % perform high cut-off
    fc(fh:fhh) = linspace(1,0,length(fh:fhh));
    
    f           = fft(dat,[],2); % FFT
    f = f.*repmat(fc,[size(f,1) 1]);
    f   = (2*real(ifft(f,[],2))); % iFFT
    filt.trial{k}=f;
end
cfg=[]; cfg.hilbert='angle'; filt=ft_preprocessing(cfg,filt);
end
%data.trial{1}=plateaufilt(data.trial{1},data.fsample,[7 13],.15);