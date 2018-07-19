
function [like_all_subj_appended, dislike_all_subj_appended] = adi_append_group_EEG (fieldtripPath)

subjList = dir(fieldtripPath);
subjList(1:2)=[]; % remove the first two ./..
    
%%   dislike:
dislike_all_subj = struct('label', [] , 'trial', [], 'time', [], 'cfg', []);

for i =  1:length(subjList) 
%     if length(dislike_all_subj) == 14
%         dislike_all_subj
%     end
    [dislike_all_subj]  =  kh_appendalldata([fieldtripPath, subjList(i,1).name, '\EEG_analysis\', '1_45Hz', '\02_interpolated\'], dislike_all_subj, 'dislike', i)
end   

size_dislike_all_subj = size(dislike_all_subj);
cfg = [];
dislike_all_subj_appended = ft_appenddata(cfg, dislike_all_subj(:,1), dislike_all_subj(:,2));

for k = 3:size_dislike_all_subj (2)
    dislike_all_subj_appended = ft_appenddata(cfg, dislike_all_subj_appended, dislike_all_subj(:,k)) ;   
end
clear dislike_all_subj

    
 %% like:

like_all_subj = struct('label', [] , 'trial', [], 'time', [], 'cfg', []);

for i =  1:length(subjList) 
    [like_all_subj]  =  kh_appendalldata( [fieldtripPath, subjList(i,1).name, '\EEG_analysis\', '1_45Hz', '\02_interpolated\'], like_all_subj, 'like', i);
end   
  
size_like_all_subj = size(like_all_subj);
cfg = [];
like_all_subj_appended = ft_appenddata(cfg, like_all_subj(:,1), like_all_subj(:,2));

for k = 3:size_like_all_subj (2)
    like_all_subj_appended = ft_appenddata(cfg, like_all_subj_appended, like_all_subj(:,k)) ;   
end

clear like_all_subj
end
    


function [condition_all_subj]=kh_appendalldata (subjFolder, condition_all_subj, name_condition, i)

files = dir(fullfile([subjFolder, name_condition,'*.mat']));
size_files = size(files);

for k = 1:(size_files(1,1))
    load ([subjFolder, files(k).name]);
     switch k 
            case 2 
               [var, cleanEEG_interp] = checkfields (var, cleanEEG_interp); 
            case 3
                [~, cleanEEG_interp] = checkfields (var(2), cleanEEG_interp); 
     end

    var(:,k) = cleanEEG_interp;
    clear cleanEEG_interp
end

cfg = [];
switch size_files(1)
    case 3
        data_allRuns = ft_appenddata(cfg, var(1), var(2), var(3))
    case 2
        data_allRuns = ft_appenddata(cfg, var(1), var(2))
end

condition_all_subj(i) = data_allRuns;

end


function [var, EEG_interp] = checkfields (var, EEG_interp)

if ~isfield(EEG_interp, 'cfg')
    EEG_interp = setfield(EEG_interp, 'cfg', []);
end

if ~isfield(var, 'cfg')
    var = setfield(var, 'cfg', []);
end

if isfield(var, 'sampleinfo')
    var = rmfield(var, 'sampleinfo');
end

if isfield(EEG_interp, 'sampleinfo')
    EEG_interp = rmfield(EEG_interp, 'sampleinfo');
end

if isfield(EEG_interp, 'dimord')
    EEG_interp = rmfield(EEG_interp, 'dimord');
end

if isfield(var, 'dimord')
    var = rmfield(var, 'dimord');
end



if isfield(EEG_interp, 'sampleinfo_orig')
    EEG_interp = rmfield(EEG_interp, 'sampleinfo_orig');
end
if isfield(var, 'sampleinfo_orig')
    var = rmfield(var, 'sampleinfo_orig');
end


% fieldnames1 = fieldnames(var);
% fieldnames2 = fieldnames(EEG_interp);
% [field, row] = setdiff(fieldnames1, fieldnames2);
% var = rmfield(var, field);
% clear field row fieldnames1 fieldnames2
fieldnames1 = fieldnames(var);
fieldnames2 = fieldnames(EEG_interp);
[field, row] = setdiff(fieldnames1, fieldnames2);
var = rmfield(var, field);
fieldnames1 = fieldnames(var);
fieldnames2 = fieldnames(EEG_interp);
[field, row] = setdiff(fieldnames2, fieldnames1);
EEG_interp = rmfield(EEG_interp, field);
EEG_interp = orderfields(EEG_interp, var);
end
