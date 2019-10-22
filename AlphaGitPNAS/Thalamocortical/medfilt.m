function filt = medfilt(data,ord)
filt = data;
for k=1:length(data.trial)
    filt.trial{k} = ft_preproc_medianfilter(data.trial{k},ord);
end