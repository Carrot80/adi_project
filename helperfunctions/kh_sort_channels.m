function [] = sort_channels(inpath, path_labels, outpath, file_name)

list_dir = dir([inpath '*.mat']);
new_label = load (path_labels);
new_label = new_label.cleanMEG.label;

for kk = 1:length(list_dir)
     
    load ([list_dir(kk).folder filesep list_dir(kk).name])
    trials = cell(1, length(cleanMEG.trial));
    chan_indx = zeros(length(new_label),1);
    for pp = 1:length(cleanMEG.label)
        chan_indx(pp) = find(strcmp(new_label{pp}, cleanMEG.label));
    end     

    for pp = 1:length(cleanMEG.trial)
        trials{pp} = cleanMEG.trial{pp}(chan_indx,:);
    end

    cleanMEG.trial = []; 
    cleanMEG.trial = trials;
    cleanMEG.label = new_label;
    save([outpath list_dir(kk).name], 'cleanMEG');        
    clear cleanMEG trials

end
















end