function [] = add_triggercode_labels(path2data, subject_list, color, triggercodes, triggercode_labels, run)


for i = 1:length(subject_list)
    
    list_clean_data = dir(fullfile([path2data subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\01_clean\alt\'], ['*500_' run '.mat'])); 
    fileList_Export_Bst2FT = dir(fullfile([path2data subject_list(i).name '\MEG_EEG_input\noisereduced\1_95Hz\02_Export_Bst2FT\'], ['*500_' run '.mat'])); 
    for k=1:length(list_clean_data)
        disp(['loading exported brainstorm file ' fileList_Export_Bst2FT(k).name ' of subject ' subject_list(i).name])
        load([path2data subject_list(i).name filesep 'MEG_EEG_input\noisereduced\1_95Hz\02_Export_Bst2FT\' fileList_Export_Bst2FT(k).name], 'RetVal')
        path_cleanMEG = [path2data subject_list(i).name filesep 'MEG_analysis\noisereduced\1_95Hz\01_clean\alt\'];
        disp(['loading cleanMEG file ' list_clean_data(k).name ' of subject ' subject_list(i).name])
        load([path_cleanMEG list_clean_data(k).name], 'cleanMEG')
        if ~isfield(cleanMEG, 'trialinfo')
            cleanMEG = setfield(cleanMEG, 'trialinfo', []);
        end
        
        if isfield (cleanMEG, 'rejectedTrials') && ~isempty(cleanMEG.rejectedTrials)
            RetVal.trial(cleanMEG.rejectedTrials) = [];
        end
        
        if ~isequal(length(cleanMEG.trial), length(RetVal.trial))
            warning('number of trials in cleanMEG and RetVal are not consistent')
        end
        for p = 1:length(cleanMEG.trial)
            fprintf('p = %d\n', p);
            if isequal(RetVal.trial{1,p}(1,:), cleanMEG.trial{1,p}(1,:)) 
                cleanMEG.trialinfo.triggerchannel(p,:) = RetVal.trial{1,p}(272,1:3052); 
                cleanMEG.trialinfo.responsechannel(p,:) = RetVal.trial{1,p}(273,1:3052); 
                cleanMEG.trialinfo.triggerlabel(p) = max(RetVal.trial{1,p}(272,950:1010));
                cleanMEG.trialinfo.response(p) = max(RetVal.trial{1,p}(273,:));
                switch cleanMEG.trialinfo.response(p)
                	case 32
                        cleanMEG.trialinfo.response_label{p} = 'dislike';
                    case 16
                          cleanMEG.trialinfo.response_label{p} = 'like';
                    case 64
                        cleanMEG.trialinfo.response_label{p} = 'dontcare';
                    otherwise
                        cleanMEG.trialinfo.response_label{p} = 'unknown';
                end
            else
               warning(['channel 1 of RetVal.trial and cleanMEG.trial in trials no. ' num2str(p) '  have different values']) 
              ind = zeros(1,length(RetVal.trial));
               for op = 1:length(RetVal.trial)
                temp(op) = isequal(cleanMEG.trial{1,p}(1,:), RetVal.trial{1,op}(1,:));
               end
               ind = find(temp);
               clear temp;
               if isempty(ind)
                  save ([path_cleanMEG 'incomplete_' list_clean_data(k).name ], 'cleanMEG' )
                  break
               end
               cleanMEG.trialinfo.triggerchannel(p,:) = RetVal.trial{1,ind}(272,:); 
               cleanMEG.trialinfo.responsechannel(p,:) = RetVal.trial{1,ind}(273,:); 
               cleanMEG.trialinfo.triggerlabel(p) = max(RetVal.trial{1,ind}(272,950:1010));
               cleanMEG.trialinfo.response(p) = max(RetVal.trial{1,ind}(273,:));
               clear ind;
               switch cleanMEG.trialinfo.response(p)
                	case 32
                        cleanMEG.trialinfo.response_label{p} = 'dislike';
                    case 16
                          cleanMEG.trialinfo.response_label{p} = 'like';
                    case 64
                        cleanMEG.trialinfo.response_label{p} = 'dontcare';
                    otherwise
                        cleanMEG.trialinfo.response_label{p} = 'unknown';
               end
            end
        end
        
        save ([path_cleanMEG list_clean_data(k).name], 'cleanMEG' )
        clear cleanMEG RetVal   
    end
        
end
    

end


function adi_recode(RetVal, cleanMEG, color, run, path_cleanMEG, condition, subject, triggercodes, triggercode_labels)

[trl_like_orig, trl_dislike_orig, trl_dontcare_orig, trl_all_orig] = adi_trialdef (subject, run);


for j=1:length(RetVal.trial)
    TriggerCodeTrl(j)= max(RetVal.trial{1,j}(272,950:1010));
end

if ~isfield(cleanMEG, 'rejectedTrials')
    cleanMEG = setfield(cleanMEG, 'rejectedTrials', []);
end

if 1 == strcmp(condition, ['dislike500_' run '.mat']) 
    try
        trl_dislike_orig(:, 5) = TriggerCodeTrl;
    catch
        if length(trl_dislike_orig)~=length(TriggerCodeTrl)
           eprimetrigger = [102 104 106 108 101 103 105 107 109]; 
           trigger_RetVal = trl_dislike_orig(:, 4); 
           ind = zeros(1, length(trigger_RetVal));
           for j=1:length(trigger_RetVal)
                [~, y] = find(trigger_RetVal(j) == eprimetrigger);
                ind(j) = str2double(char(triggercodes{y}));
           end
           trl_dislike_orig(:, 5) = ind;
           TriggerCodeTrl = TriggerCodeTrl';
           trl_dislike_orig(4,:)=[];
           if 1==~any(1==ismember(TriggerCodeTrl, trl_dislike_orig(:,5)))
               error('bitte trigger codes überprüfen')
           end
        end
    end

    trl_dislike_badTrials = trl_dislike_orig(cleanMEG.rejectedTrials,:);
    cleanMEG.trialinfo.trl_dislike_orig = trl_dislike_orig;
    cleanMEG.trialinfo.trl_dislike_badTrials = trl_dislike_badTrials;
    cleanMEG.trialinfo.rejected_trials = cleanMEG.rejectedTrials;
    trl_dislike = trl_dislike_orig;
    trl_dislike(cleanMEG.rejectedTrials,:) = [];
    cleanMEG.trialinfo.sampleinfo_inclusive_new_dontcare = trl_dislike;
    
else
    try
    trl_like_orig(:, 5) = TriggerCodeTrl;
    catch
        warning(['Triggercodes überprüfen'])
        if length(trl_like_orig)~=length(TriggerCodeTrl)
           eprimetrigger = [102 104 106 108 101 103 105 107 109]; 
           trigger_RetVal = trl_like_orig(:, 4); 
           ind = zeros(1, length(trigger_RetVal));
           for j=1:length(trigger_RetVal)
                [~, y] = find(trigger_RetVal(j) == eprimetrigger);
                ind(j) = str2double(char(triggercodes{y}));
           end
           trl_like_orig(:, 5) = ind;
           TriggerCodeTrl = TriggerCodeTrl';
           TriggerCodeTrl(4)=[];
           if ~all(1==ismember(TriggerCodeTrl, trl_like_orig(:,5)))
               error('bitte trigger codes überprüfen')
           end
        end
    end
    trl_like_badTrials = trl_like_orig(cleanMEG.rejectedTrials,:);
    cleanMEG.trialinfo.trl_like_orig = trl_like_orig;
    cleanMEG.trialinfo.trl_like_badTrials = trl_like_badTrials;
    cleanMEG.trialinfo.rejected_trials = cleanMEG.rejectedTrials;
    trl_like = trl_like_orig;
    trl_like(cleanMEG.rejectedTrials,:) = [];
    cleanMEG.trialinfo.sampleinfo_inclusive_new_dontcare = trl_like;
end
    
TriggerCodeTrl(cleanMEG.rejectedTrials) = [];
[~, indTrl] = find(ismember(TriggerCodeTrl, str2double(triggercodes2recode)));

if 1 == strcmp(condition, ['dislike500_' run '.mat']) 
    cleanMEG.trialinfo.sampleinfo_new_dontcare = trl_dislike(indTrl, :);
    cleanMEG.trialinfo.sampleinfo = trl_dislike;
else
    cleanMEG.trialinfo.sampleinfo_new_dontcare = trl_like(indTrl, :);
    cleanMEG.trialinfo.sampleinfo = trl_like;
end

cleanMEG.trialinfo.sampleinfo(indTrl, :) = [];

% integrate trials into dontcare variable:
if exist([path_cleanMEG 'dontcare500_' run '.mat'], 'file')
    dontcare_integ = load ([path_cleanMEG 'dontcare500_' run '.mat']);
    if any(ismember(cleanMEG.trialinfo.sampleinfo_new_dontcare(:,1), dontcare_integ.cleanMEG.trialinfo.sampleinfo(:,1)))
        return
    end
elseif exist([path_cleanMEG 'alt\dontcare500_' run '.mat'], 'file')  
    dontcare_integ = load ([path_cleanMEG 'alt\dontcare500_' run '.mat']);  
else 
    dontcare_integ.cleanMEG = struct;
    dontcare_integ.cleanMEG.trial = {};
    dontcare_integ.cleanMEG.time = {};
    dontcare_integ.cleanMEG.sampleinfo = [];
    dontcare_integ.cleanMEG.ChannelFlag_Bst = {};
    dontcare_integ.cleanMEG.grad = cleanMEG.grad;
end

if ~isfield(dontcare_integ.cleanMEG, 'rejectedTrials')
    dontcare_integ.cleanMEG = setfield(dontcare_integ.cleanMEG, 'rejectedTrials', []);
end

dontcare_integ.cleanMEG.trial = [dontcare_integ.cleanMEG.trial cleanMEG.trial(indTrl)]; 
dontcare_integ.cleanMEG.time = [dontcare_integ.cleanMEG.time cleanMEG.time(indTrl)]; 
dontcare_integ.cleanMEG.sampleinfo = [dontcare_integ.cleanMEG.sampleinfo; cleanMEG.sampleinfo(indTrl,:)]; 
dontcare_integ.cleanMEG.ChannelFlag_Bst = [dontcare_integ.cleanMEG.ChannelFlag_Bst cleanMEG.ChannelFlag_Bst(indTrl)]; 
dontcare_integ.cleanMEG.trialinfo.trl_dontcare_orig = trl_dontcare_orig;
dontcare_integ.cleanMEG.trialinfo.trl_dontcare_badTrials = trl_dontcare_orig(dontcare_integ.cleanMEG.rejectedTrials, :);

trl_dontcare = trl_dontcare_orig;
trl_dontcare(:,5) = 0;
trl_dontcare(dontcare_integ.cleanMEG.rejectedTrials,:) = [];

if 1 == strcmp(condition, ['dislike500_' run '.mat']) 
    if ~isfield(dontcare_integ.cleanMEG.trialinfo, 'sampleinfo')
        dontcare_integ.cleanMEG.trialinfo.sampleinfo = [trl_dontcare; trl_dislike(indTrl,:)];
    else
        if ~any(ismember(trl_dislike(indTrl,1), dontcare_integ.cleanMEG.trialinfo.sampleinfo(:,1))) 
            dontcare_integ.cleanMEG.trialinfo.sampleinfo = [dontcare_integ.cleanMEG.trialinfo.sampleinfo; trl_dislike(indTrl,:)];
        end
    end
else
    if ~isfield(dontcare_integ.cleanMEG.trialinfo, 'sampleinfo')
        dontcare_integ.cleanMEG.trialinfo.sampleinfo = [trl_dontcare; trl_like(indTrl,:)];
    else
        if ~any(ismember(trl_like(indTrl,1), dontcare_integ.cleanMEG.trialinfo.sampleinfo(:,1) ))
            dontcare_integ.cleanMEG.trialinfo.sampleinfo = [dontcare_integ.cleanMEG.trialinfo.sampleinfo; trl_like(indTrl,:)];
        end
    end
end  


%%

% clean cleanMEG from new dontcare trials
cleanMEG.trial(indTrl) = [];
for j=1:length(indTrl)
    cleanMEG.time{1, indTrl(j)} = [];
end
cleanMEG.time(indTrl)=[];
cleanMEG.ChannelFlag_Bst(indTrl)=[];
cleanMEG.sampleinfo(indTrl,:)=[];

% save new variables:
save ([path_cleanMEG condition], 'cleanMEG')
clear cleanMEG
cleanMEG = dontcare_integ.cleanMEG;
save ([path_cleanMEG 'dontcare500_' run '.mat'], 'cleanMEG' )

end
