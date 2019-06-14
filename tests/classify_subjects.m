
subj1 = 'nl_adi_06';
subj2 = 'nl_adi_08';
path_subj1 = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\' subj1 '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
path_subj2 = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\' subj2 '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];

files_subj1 = dir(fullfile(path_subj1, '*.mat'));

for p = 1:length(files_subj1)
    load([files_subj1(p).folder filesep files_subj1(p).name])

    while p > 1

        [var, cleanMEG_interp] = checkfields (var(p-1), cleanMEG_interp); 

    end


    for k=1:length(cleanMEG_interp.trial)
        cleanMEG_interp.run(k) = str2double(files_subj1(p).name(end-4));
        cleanMEG_interp.subj(k) = str2double(subj1(end-1:end));
        cleanMEG_interp.subj_label(k) = 1;
        if 1 == contains(files_subj1(p).name, 'Dislike')
            cleanMEG_interp.condition{k} = 'dislike';
            cleanMEG_interp.cond_label(k) = 2;
        elseif 1 == contains(files_subj1(p).name, 'Neu_Like')
            cleanMEG_interp.condition{k} = 'like';
            cleanMEG_interp.condition(k) = 1;
        end

    end

    var(:,p) = cleanMEG_interp;
    clear cleanMEG_interp

end

%%
function [var, MEG_interp] = checkfields (var, MEG_interp)

    if ~isfield(MEG_interp, 'dimord')
        MEG_interp = setfield(MEG_interp, 'dimord', 'chan_time');
    end
    if ~isfield(var, 'dimord')
        var = setfield(var, 'dimord', 'chan_time');
    end
    if ~isfield(MEG_interp, 'sampleinfo')
        MEG_interp.sampleinfo = MEG_interp.sampleinfo_orig;
    end
    if ~isfield(var, 'sampleinfo')
        var.sampleinfo = var.sampleinfo_orig;
    end

    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEG_interp);
    [field, row] = setdiff(fieldnames1, fieldnames2);
    var = rmfield(var, field);
    clear field row fieldnames1 fieldnames2
    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEG_interp);
    [field, row] = setdiff(fieldnames2, fieldnames1);
    MEG_interp = rmfield(MEG_interp, field);
    MEG_interp = orderfields(MEG_interp, var);
    
end