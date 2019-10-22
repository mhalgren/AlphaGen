clear
%subset of subjects with HGP
miall=[];
for k=[1 2 3 10 13 14 15 51 52]
    load(strcat('L',num2str(k),'Resultsv10.mat'),'lbl','mod_ind','SR')
    if length(SR)==length(diag(mod_ind))
        if isempty(miall)
        miall=diag(mod_ind)'; srall=SR; 
        else
            miall=[miall diag(mod_ind)']; srall=[srall SR];
        end
    end
end

[r,pv]=corr(miall',srall')