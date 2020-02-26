function [] = count_runs(mainpath, subjects_dir, filename)

for ii = 1:length(subjects_dir)
    
   dir_data = dir([mainpath subjects_dir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\' '*.mat']);
   for kk = 1:length(dir_data)
        load ([dir_data(kk).folder filesep dir_data(kk).name], 'cleanMEG_interp')
        for pp = 1:length(cleanMEG_interp.trial)
            cleanMEG_interp.trialinfo.run{pp} = dir_data(kk).name(end-4);
        end

        [data_bpfreq_res_sel] = adi_bpfilter(cleanMEG_interp, 'bp1_45Hz');

        data(kk) = cleanMEG_interp;
        clear cleanMEG_interp data_bpfreq_res_sel
    end

     
    
end







end