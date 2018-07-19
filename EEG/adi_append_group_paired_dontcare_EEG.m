
function [group_data_like, group_data_dislike, group_data_dontcare] = adi_append_group_paired_dontcare_EEG (fieldtripPath)
    
subjList = dir( fieldtripPath );
subjList(1:2)=[]; % remove the first two ./..

num_dislike = 1;
num_like = 1;
num_dontcare = 1;

dislike_all_subj = struct('label', [] ,  'dimord', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', []);
like_all_subj = struct('label', [] ,  'dimord', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', [] );
dontcare_all_subj = struct('label', [] ,  'dimord', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', []);

for i =  1:length(subjList) 
    inPath = ([fieldtripPath, subjList(i,1).name, '\EEG_analysis\1_45Hz\02_interpolated\']);
    [num_dislike, num_like, num_dontcare, dislike_all_subj, like_all_subj, dontcare_all_subj]  =  kh_appendalldata(inPath, num_dislike, num_like, num_dontcare, dislike_all_subj, like_all_subj, dontcare_all_subj);
end   
    
%% dislike:
size_dislike_all_subj = size(dislike_all_subj);
cfg = [];
dislike_appended_paired = ft_appenddata(cfg, dislike_all_subj(:,1), dislike_all_subj(:,2));

for k = 3:size_dislike_all_subj (2)
    dislike_appended_paired = ft_appenddata(cfg, dislike_appended_paired, dislike_all_subj(:,k)) ;   
end

    
    %% like:
 
size_like_all_subj=size(like_all_subj);
cfg=[];
like_appended_paired = ft_appenddata(cfg, like_all_subj(:,1), like_all_subj(:,2));

for k = 3:size_like_all_subj (2)
    like_appended_paired = ft_appenddata(cfg, like_appended_paired, like_all_subj(:,k))  ;  
end

    %% dontcare:

size_dontcare_all_Runs = size(dontcare_all_subj);
cfg = [];
dontcare_appended_paired = ft_appenddata(cfg, dontcare_all_subj(1,:), dontcare_all_subj(2,:));

for k = 3:size_dontcare_all_Runs (1)
    dontcare_appended_paired = ft_appenddata(cfg, dontcare_appended_paired, dontcare_all_subj(k,:)) ;   
end
    

end



function [num_dislike, num_like, num_dontcare, dislike_all_subj, like_all_subj, dontcare_all_subj] = kh_appendalldata ( inPath, num_dislike, num_like, num_dontcare, dislike_all_subj, like_all_subj, dontcare_all_subj)

%% dontcare:

files_dontcare = dir(fullfile([inPath, '*dontcare*.mat']));
size_files_dontcare = size(files_dontcare);

for k = 1:(size_files_dontcare(1,1))
    load ([subjFolder, files_dontcare(k).name]);
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
switch size_files_dontcare(1)
    case 3
        dontcare_allRuns = ft_appenddata(cfg, var(1), var(2), var(3))
    case 2
        dontcare_allRuns = ft_appenddata(cfg, var(1), var(2))
end


fieldnames_dontcare = fieldnames(dontcare_allRuns);
fieldnames_dontcare_all_subj = fieldnames(dontcare_all_subj);
diff_fieldnames = setdiff (fieldnames_dontcare, fieldnames_dontcare_all_subj);
if ~isempty(diff_fieldnames)
    dontcare_allRuns = rmfield(dontcare_allRuns, diff_fieldnames);
end
dontcare_allRuns = orderfields(dontcare_allRuns, dontcare_all_subj);
dontcare_all_subj(num_dontcare,:) = dontcare_allRuns;
nr_trl_dontcare_allRuns = length(dontcare_allRuns.trial);
clear dontcare_allRuns     

num_dontcare = num_dontcare+1;
   

%%
    file_dislike = ([path_appended, 'dislike_allRuns_' freqname, '.mat']);
    if exist (file_dislike, 'file')
       load (file_dislike);
       if length(dislike_allRuns.trial) > nr_trl_dontcare_allRuns
            dislike_allRuns.trial = dislike_allRuns.trial(1:nr_trl_dontcare_allRuns);
            dislike_allRuns.time = dislike_allRuns.time(1:nr_trl_dontcare_allRuns);
             if isfield (dislike_allRuns, 'sampleinfo')
                dislike_allRuns.sampleinfo = dislike_allRuns.sampleinfo(1:nr_trl_dontcare_allRuns,:);
             else
               dislike_allRuns = setfield(dislike_allRuns, 'sampleinfo', []);
             end
       else   
            if ~isfield (dislike_allRuns, 'sampleinfo')
              dislike_allRuns = setfield(dislike_allRuns, 'sampleinfo', []);
            end
       end
       if ~isfield (dislike_allRuns, 'dimord')
           dislike_allRuns = setfield(dislike_allRuns, 'dimord', []);
       end
       fieldnames_dislike = fieldnames(dislike_allRuns);
       fieldnames_dislike_all_subj = fieldnames(dislike_all_subj);
       diff_fieldnames = setdiff (fieldnames_dislike, fieldnames_dislike_all_subj);
       if ~isempty(diff_fieldnames)
            dislike_allRuns = rmfield(dislike_allRuns, diff_fieldnames);
       end
       dislike_allRuns = orderfields(dislike_allRuns, dislike_all_subj);    
       if num_dislike == 1
            dislike_all_subj = dislike_allRuns;
       else
            dislike_all_subj(num_dislike)=dislike_allRuns;
       end
        
       clear dislike_allRuns        

   
        num_dislike = num_dislike+1;
    end
 %% like:   
  
  file_like = ([path_appended, 'like_allRuns_', freqname, '.mat']);
  if exist (file_like, 'file')
    load(file_like);
    if length(like_allRuns.trial) > nr_trl_dontcare_allRuns    
        like_allRuns.trial = like_allRuns.trial(1:nr_trl_dontcare_allRuns);
        like_allRuns.time = like_allRuns.time(1:nr_trl_dontcare_allRuns);
        if isfield(like_allRuns, 'sampleinfo')
             like_allRuns.sampleinfo = like_allRuns.sampleinfo(1:nr_trl_dontcare_allRuns,:);
        else
           like_allRuns = setfield(like_allRuns, 'sampleinfo', []);
        end
        if ~isfield(like_allRuns, 'dimord')
             like_allRuns = setfield(like_allRuns, 'dimord', []);
        end
    else
        if ~isfield (like_allRuns, 'sampleinfo')
              like_allRuns = setfield(like_allRuns, 'sampleinfo', []);
        end
        if ~isfield(like_allRuns, 'dimord')
             like_allRuns = setfield(like_allRuns, 'dimord', []);
        end
    end

    fieldnames_like = fieldnames(like_allRuns);
    fieldnames_like_all_subj = fieldnames(like_all_subj);
    diff_fieldnames = setdiff (fieldnames_like, fieldnames_like_all_subj);
       if ~isempty(diff_fieldnames)
            like_allRuns = rmfield(like_allRuns, diff_fieldnames);
       end
    like_allRuns = orderfields(like_allRuns, like_all_subj);    
    if num_like == 1
        like_all_subj = like_allRuns;
    else
        like_all_subj(num_like) = like_allRuns;
    end

    clear like_allRuns     
 
    num_like = num_like+1;
  end
end



function [condition_all_subj] = kh_appendalldata (subjFolder, condition_all_subj, name_condition, i)

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

if ~isfield(EEG_interp, 'dimord')
    EEG_interp = setfield(EEG_interp, 'dimord', 'chan_time');
end
if ~isfield(var, 'dimord')
    var = setfield(var, 'dimord', 'chan_time');
end
if ~isfield(EEG_interp, 'sampleinfo')
    EEG_interp.sampleinfo = EEG_interp.sampleinfo_orig;
end
if ~isfield(var, 'sampleinfo')
    var.sampleinfo = var.sampleinfo_orig;
end

fieldnames1 = fieldnames(var);
fieldnames2 = fieldnames(EEG_interp);
[field, row] = setdiff(fieldnames1, fieldnames2);
var = rmfield(var, field);
clear field row fieldnames1 fieldnames2
fieldnames1 = fieldnames(var);
fieldnames2 = fieldnames(EEG_interp);
[field, row] = setdiff(fieldnames2, fieldnames1);
EEG_interp = rmfield(EEG_interp, field);
EEG_interp = orderfields(EEG_interp, var);
end

