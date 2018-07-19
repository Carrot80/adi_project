function [data_bpfreq] = adi_bpfilter_EEG(inPath, outPath, bpname)


switch bpname
    case 'delta'
        bpfreq = 4;
    case 'theta'
        bpfreq = [4 8];
    case 'alpha'
        bpfreq = [8 13];
    case 'beta'
        bpfreq = [13 25];
    case 'low_gamma'
        bpfreq = [25 45];
end

cfg=[];
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') || 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

list = dir(fullfile([inPath, '*.mat'])); 
for k=1:length(list)
    if ~exist([outPath, list(k).name(1:end-4), '_', bpname, '.mat'], 'file')
        load([inPath, list(k).name]); 
        [data_bpfreq] = ft_preprocessing(cfg, cleanEEG_interp);
        save ([outPath, list(k).name(1:end-4), '_', bpname, '.mat'], 'data_bpfreq')
        clear data_bpfreq cleanEEG_interp
    end
end
     
