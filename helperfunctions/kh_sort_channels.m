function [] = sort_channels(inpath, path_labels, outpath, file_name)

% noch anpassen, funktioniert noch nicht
list_dir = dir([inpath '*.mat']);
new_label = load (path_labels);
label = cleanMEG.label;
clear cleanMEG

for kk = 1:length(list_dir)
     
    load ([list_dir(kk).folder filesep list_dir(kk).name])
    trials = cell(1, length(cleanMEG.trial));
    chan_indx = zeros(length(label),1);
    for pp = 1:length(cleanMEG.label)
        chan_indx(pp) = find(strcmp(label{pp}, cleanMEG.label));
    end     

    for pp = 1:length(cleanMEG.trial)
        trials{pp} = cleanMEG.trial{pp}(chan_indx,:);
    end

    cleanMEG.trial = []; 
    cleanMEG.trial = trials;
    cleanMEG.label = label;
    save(['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\nl_adi_04\MEG_analysis\noisereduced\1_95Hz\01_clean\' list_dir(kk).name], 'cleanMEG');        
    clear cleanMEG trials

end
















end