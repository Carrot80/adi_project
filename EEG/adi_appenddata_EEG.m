
function  [dislike_allRuns, like_allRuns] = adi_appenddata_EEG(inPath, filter)

% append Runs 

% dislike
files_dislike = dir(fullfile(inPath, 'dislike*.mat'));
size_files_dislike = size(files_dislike);

for i = 1:(size_files_dislike(1,1))
    load ([inPath, files_dislike(i).name]);
     switch i 
            case 2 
               [var, cleanEEG_interp] = checkfields (var, cleanEEG_interp); 
            case 3
                [~, cleanEEG_interp] = checkfields (var(2), cleanEEG_interp); 
     end
    var(:,i) = cleanEEG_interp;
    clear cleanEEG_interp
end

cfg = [];
switch size_files_dislike(1)
    case 3
        dislike_allRuns = ft_appenddata(cfg, var(1), var(2), var(3))
    case 2
        dislike_allRuns = ft_appenddata(cfg, var(1), var(2))
end
clear var    

%% like 

files_like = dir(fullfile(inPath, 'like*.mat'));
size_files_like = size(files_like);

for i = 1:(size_files_like(1,1))
    load ([inPath, files_like(i).name]);
    switch i 
        case 2 
           [var, cleanEEG_interp] = checkfields (var, cleanEEG_interp); 
        case 3
            [~, cleanEEG_interp] = checkfields (var(2), cleanEEG_interp); 
    end
    var(:,i)= cleanEEG_interp;
    clear cleanEEG_interp;           
end

cfg = [];
switch size_files_like(1) 
    case 3
        like_allRuns = ft_appenddata(cfg, var(1), var(2), var(3))
    case 2
        like_allRuns = ft_appenddata(cfg, var(1), var(2))
end

clear var 
 
%% dontcare
%     if ~exist (strcat(outPath, 'dontcare_allRuns_', filter, '.mat'))
%         files_dontcare = dir(fullfile(inPath, 'dont*.mat'));
%         size_files_dontcare = size(files_dontcare);
% 
%         if size_files_dontcare(1,1) == 0
%             return
%         else
%             for i=1:(size_files_dontcare(1,1))
%                 load (strcat(inPath, files_dontcare(i).name))
%             switch i 
%                 case 2 
%                    [var, cleanEEG_interp] = checkfields (var, cleanEEG_interp); 
%                 case 3
%                     [~, cleanEEG_interp] = checkfields (var(2), cleanEEG_interp); 
%             end
%                    var(:,i)= cleanEEG_interp;
%                    clear cleanEEG_interp;      
%             end
% 
%             cfg=[];
%             switch size_files_dontcare(1)
%                 case 3
%                     dontcare_allRuns = ft_appenddata(cfg, var(1), var(2), var(3))
%                     save (strcat(outPath, filesep, 'dontcare_allRuns_', filter, '.mat'), 'dontcare_allRuns');
% 
%                 case 2
%                     dontcare_allRuns = ft_appenddata(cfg, var(1), var(2))
%                     save (strcat(outPath, filesep, 'dontcare_allRuns_', filter, '.mat'), 'dontcare_allRuns');
%                 case 1
%                     dontcare_allRuns = var; 
%                     save (strcat(outPath, 'dontcare_allRuns_', filter, '.mat'), 'dontcare_allRuns');
%             end
%         end
%     end
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