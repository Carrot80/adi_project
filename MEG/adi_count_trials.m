function [] = adi_count_trials(subjlist)


for ii = 1:length(subjlist)
    
    filelist = dir([subjlist(ii).folder filesep, subjlist(ii).folder filesep 'nl_adi_04\MEG_analysis\noisereduced\1_95Hz\01_clean\']);
    for kk = 1:length(filelist)
        
        load ([filelist(kk).folder filesep filelist(kk).name])
        count.like
        
        
        % hier weitermachen
        
        
        
        
    end
















end