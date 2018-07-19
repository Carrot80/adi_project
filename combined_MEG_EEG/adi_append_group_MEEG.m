
function  [group_data_like, group_data_dislike] = adi_appenddata_MEEG(inPathEEG, EEG, inPathMEG, MEG, subj, group_data_like, group_data_dislike, i)

% append Runs for MEG and  EEG

[EEG_allRuns_dislike] = load_MEEG(inPathEEG, EEG, 'dislike');
[MEG_allRuns_dislike] = load_MEEG(inPathMEG, MEG, 'dislike');
[EEG_allRuns_like] = load_MEEG(inPathEEG, EEG, 'like');
[MEG_allRuns_like] = load_MEEG(inPathMEG, MEG, 'like');

[data_combined_MEG_EEG_dislike] = append_data_MEEG(MEG_allRuns_dislike, EEG_allRuns_dislike, subj);
[data_combined_MEG_EEG_like] = append_data_MEEG(MEG_allRuns_like, EEG_allRuns_like, subj);

group_data_like(i).label = data_combined_MEG_EEG_like.label;
group_data_like(i).trial = data_combined_MEG_EEG_like.trial;
group_data_like(i).time = data_combined_MEG_EEG_like.time;
group_data_like(i).cfg = data_combined_MEG_EEG_like.cfg;

group_data_dislike(i).label = data_combined_MEG_EEG_dislike.label;
group_data_dislike(i).trial = data_combined_MEG_EEG_dislike.trial;
group_data_dislike(i).time = data_combined_MEG_EEG_dislike.time;
group_data_dislike(i).cfg = data_combined_MEG_EEG_dislike.cfg;

end

function [var] = load_MEEG(inPath, MEEG, condition)


    files = dir(fullfile([inPath, condition, '*.mat']));
    size_files = size(files);

    for i = 1:(size_files(1,1))
        load ([inPath, files(i).name]);
        if 1 == strcmp(MEEG, 'MEG')
            cleanMEEG_interp = cleanMEG_interp;
            clear cleanMEG_interp
        elseif 1 == strcmp(MEEG, 'EEG')
            cleanMEEG_interp = cleanEEG_interp;
            clear cleanEEG_interp
        end
         switch i 
                case 2 
                   [var, cleanMEEG_interp] = checkfields (var, cleanMEEG_interp, 'i_eq_2'); 
                case 3
                    [~, cleanMEEG_interp] = checkfields (var(2), cleanMEEG_interp, 'i_eq_3'); 
         end
        var(:,i) = cleanMEEG_interp;
        clear cleanMEEG_interp
    end


end

function [combined_MEG_EEG_data] = append_data_MEEG(MEG_data, EEG_data, subj)

if 1 == strcmp (subj, 'nl_adi_07')
    MEG_data(2) = [];
end

for k = 1:length(MEG_data)
    
 % für MEG sampleinfo aufbauen, sofern keine vorhanden ist: 

    if ~isfield(MEG_data(k), 'sampleinfo')
         MEG_data(k).sampleinfo = MEG_data(k).sampleinfo_orig;
        if length(MEG_data(k).trial) ~= length(MEG_data(k).sampleinfo_orig) 
            MEG_data(k).sampleinfo(MEG_data(k).rejectedTrials, :) = [];
        end
    end
    
     % um MEG und EEG zusammenzufügen, müssen gleiche Anzahl an Trials vorhanden
    % sein:
   
    rejected_trials_in_EEG_2reject_in_MEG = setdiff(MEG_data(k).sampleinfo(:,1), EEG_data(k).sampleinfo(:,1));
    rejected_trials_in_MEG_2reject_in_EEG = setdiff(EEG_data(k).sampleinfo(:,1), MEG_data(k).sampleinfo(:,1));
    
    if length(rejected_trials_in_EEG_2reject_in_MEG) > 10 
%         [sampleinfo] = mk_sampleinfo(MEG_data)
    end
    if length(rejected_trials_in_MEG_2reject_in_EEG) > 10 
        [sampleinfo] = mk_sampleinfo(MEG_data)
    end
    
    if ~isempty(rejected_trials_in_EEG_2reject_in_MEG)
        for m = 1:length(rejected_trials_in_EEG_2reject_in_MEG)
            [ind(m)] = find (MEG_data(k).sampleinfo(:,1) == rejected_trials_in_EEG_2reject_in_MEG(m))
        end
        for m = 1:length(rejected_trials_in_EEG_2reject_in_MEG)
            MEG_data(k).trial{1, ind(m)} = [];
            MEG_data(k).time{1, ind(m)} = [];
        end
        MEG_data(k).trial(cellfun('isempty',MEG_data(k).trial)) = [];
        MEG_data(k).time(cellfun('isempty',MEG_data(k).time)) = [];
        MEG_data(k).sampleinfo(ind',:) = [];
        clear ind rejected_trials_in_EEG_2reject_in_MEG
    end
    if ~isempty(rejected_trials_in_MEG_2reject_in_EEG)
         for m = 1:length(rejected_trials_in_MEG_2reject_in_EEG)
            [ind(m)] = find (EEG_data(k).sampleinfo(:,1) == rejected_trials_in_MEG_2reject_in_EEG(m))
         end
         for m = 1:length(rejected_trials_in_MEG_2reject_in_EEG)
            EEG_data(k).trial{1, ind(m)} = [];
            EEG_data(k).time{1, ind(m)} = [];
         end
         EEG_data(k).trial(cellfun('isempty',EEG_data(k).trial)) = [];
         EEG_data(k).time(cellfun('isempty',EEG_data(k).time)) = [];
         clear ind rejected_trials_in_MEG_2reject_in_EEG
    end 
end


    % Referenzsensoren rauswerfen:
    if isfield (MEG_data,  'sampleinfo')
        MEG_data = rmfield(MEG_data, 'sampleinfo')
    end
    if isfield (MEG_data,  'sampleinfo_orig')
        MEG_data = rmfield(MEG_data, 'sampleinfo_orig')
    end
    if isfield(EEG_data, 'sampleinfo')
    	EEG_data = rmfield(EEG_data, 'sampleinfo')
    end
    if isfield(EEG_data, 'sampleinfo_orig')
        EEG_data = rmfield(EEG_data, 'sampleinfo_orig')
    end
    
    if isfield(MEG_data, 'ChannelFlag_Bst')
        MEG_data = rmfield(MEG_data, 'ChannelFlag_Bst')
    end
    cfg = [];
    cfg.channel = 'meg';
    MEG_data_red = ft_selectdata(cfg, MEG_data(1))
    for p = 2:length(MEG_data)
        MEG_data_red(p) = ft_selectdata(cfg, MEG_data(p))
    end
    
    clear MEG_data
    
    cfg = [];
    switch length(MEG_data_red)
        case 3
            MEG_allRuns_appended = ft_appenddata(cfg, MEG_data_red(1), MEG_data_red(2), MEG_data_red(3))
        case 2
            MEG_allRuns_appended = ft_appenddata(cfg, MEG_data_red(1), MEG_data_red(2))
    end
    clear MEG_data_red
    
    % MEG muss von 1:45Hz gefilert werden wie EEG
    cfg = [];
    cfg.trials        = 'all'; 
    cfg.feedback      = 'yes';
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = 45;
    [MEG_allRuns_appended_bpfreq] = ft_preprocessing(cfg, MEG_allRuns_appended);   
    clear MEG_allRuns_appended    
    
    if isfield(MEG_allRuns_appended_bpfreq, 'dimord')
        MEG_allRuns_appended_bpfreq = rmfield(MEG_allRuns_appended_bpfreq, 'dimord')
    end
 
    cfg = [];
    switch length(EEG_data)
        case 3
            EEG_allRuns_appended = ft_appenddata(cfg, EEG_data(1), EEG_data(2), EEG_data(3))
        case 2
            EEG_allRuns_appended = ft_appenddata(cfg, EEG_data(1), EEG_data(2))
    end
    clear EEG_data
    
    if ~isfield(EEG_allRuns_appended, 'fsample')
        EEG_allRuns_appended = setfield(EEG_allRuns_appended, 'fsample', 1.017250000000000e+03)
    end
    
    fieldnames_MEG = fieldnames(MEG_allRuns_appended_bpfreq);
    fieldnames_EEG = fieldnames(EEG_allRuns_appended);
    diff = setdiff(fieldnames_MEG, fieldnames_EEG);
    MEG_allRuns_appended_bpfreq = rmfield(MEG_allRuns_appended_bpfreq, diff);
    
    fieldnames_MEG = fieldnames(MEG_allRuns_appended_bpfreq);
    fieldnames_EEG = fieldnames(EEG_allRuns_appended);
    diff = setdiff(fieldnames_EEG, fieldnames_MEG);
    EEG_allRuns_appended = rmfield(EEG_allRuns_appended, diff);
    MEG_allRuns_appended_bpfreq = orderfields(MEG_allRuns_appended_bpfreq, EEG_allRuns_appended)
    
    if ~isequal(MEG_allRuns_appended_bpfreq.fsample, EEG_allRuns_appended.fsample)
        MEG_allRuns_appended_bpfreq.fsample =  EEG_allRuns_appended.fsample;
    end
    

    for n = 1:length(MEG_allRuns_appended_bpfreq.trial) 
        trial_samples(n) = length(MEG_allRuns_appended_bpfreq.trial{1,n})
    end
    
    ind_samples = find( trial_samples == 3053);
    ind_samples_good = find( trial_samples ==  3052);
    
    if ~isempty (ind_samples)
        for n = 1:length(ind_samples)
            MEG_allRuns_appended_bpfreq.time{1,ind_samples(n)} = MEG_allRuns_appended_bpfreq.time{ind_samples_good(1)};
            MEG_allRuns_appended_bpfreq.trial{1,ind_samples(n)}(:,1) = [];
        end
    end
        
    % Fehlermeldung kommt im folgenden, da data.time ab der
    % 15.Nachkommastelle nicht übereinstimmt, deshalb data.time anpassen:
    
   
  
    EEG_allRuns_appended.time = MEG_allRuns_appended_bpfreq.time;
    cfg = [];
    cfg.appenddim = 'chan';
    combined_MEG_EEG_data = ft_appenddata(cfg, MEG_allRuns_appended_bpfreq, EEG_allRuns_appended)
    
end

function [var, MEEG_interp] = checkfields (var, MEEG_interp, i_eq)

if isfield(var, 'fsample')
   var = rmfield(var, 'fsample');
end

if isfield(MEEG_interp, 'fsample')
   MEEG_interp = rmfield(MEEG_interp, 'fsample');
end

if 1 == strcmp(i_eq, 'i_eq_2') 
    if ~isfield(var, 'dimord')
        var = setfield(var, 'dimord', 'chan_time');
    end
    if ~isfield(var, 'sampleinfo')
         if isfield(var, 'sampleinfo_orig')
            var.sampleinfo = var.sampleinfo_orig;
         else 
            [sampleinfo] = mk_sampleinfo(var)
             var.sampleinfo = sampleinfo;
         end
    end
    
end
        if ~isfield(MEEG_interp, 'dimord')
            MEEG_interp = setfield(MEEG_interp, 'dimord', 'chan_time');
        end
        if ~isfield(MEEG_interp, 'sampleinfo')
            if isfield(MEEG_interp, 'sampleinfo_orig')
                MEEG_interp.sampleinfo = MEEG_interp.sampleinfo_orig;
            else 
                 [sampleinfo] = mk_sampleinfo(MEEG_interp)
                MEEG_interp.sampleinfo = sampleinfo;
            end
        end

    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEEG_interp);
    [field, row] = setdiff(fieldnames1, fieldnames2);
   
  
    if 1 == strcmp(i_eq, 'i_eq_2')    
        var = rmfield(var, field);
    else 
        for p = 1:length(field)
            MEEG_interp = setfield(MEEG_interp, field{p} , '[]')
        end
    end
    
    
    clear field row fieldnames1 fieldnames2
    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEEG_interp);
    [field, row] = setdiff(fieldnames2, fieldnames1);
    MEEG_interp = rmfield(MEEG_interp, field);
    MEEG_interp = orderfields(MEEG_interp, var);
    
end


function [sampleinfo] = mk_sampleinfo(data)

len = length(data.trial)
sampleinfo = zeros(len,2);
sampleinfo(1, 1) = 1;
sampleinfo(1, 2) = length(data.trial{1,1});
o = 1;
for p = 2:len
    sampleinfo(p,1) = sampleinfo(1,1)+sampleinfo(o,2) % 64 statt 65
    sampleinfo(p,2) = sampleinfo(p,1)+sampleinfo(1,2)-1
    o = o+1;
end

end