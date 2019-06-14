function [deviance] = mean_IOD (path2data, subject, deviance, subjcounter)

files = dir([path2data '*.mat']);

for k=1:length(files)
    load ([files(k).folder filesep files(k).name])
    
    for p=1:size(cleanMEG.trialinfo.triggerchannel,1)
        ind_trig = find(cleanMEG.trialinfo.triggerchannel(p,:));
        ind_IOD = find(cleanMEG.trialinfo.triggerchannel(p,:)==512);
        sample_first_trig = ind_trig(1);
        sample_first_IOD = ind_IOD(1);
        diff=sample_first_IOD-sample_first_trig;
        deviance(subjcounter).file(k).diff(p) = diff;
        clear diff 
    end
    
    clear cleanMEG
    
end


















end