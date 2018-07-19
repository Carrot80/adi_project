
function [like_all_subj_appended, dislike_all_subj_appended] = adi_appendgroup_MEG_sourcespace (path2disc, freq, cfg_virtsens)

subjList = dir(path2disc);
subjList(1:2)=[]; % remove the first two ./..
    
%%   dislike:
dislike_all_subj = struct('trial', [], 'time', []);

for i =  1:length(subjList) 
    load ([path2disc subjList(i).name '\runs_appended\virtsens\' cfg_virtsens '_', 'dislike_allRuns_' freq '.mat'], 'vs_allRuns');
    dislike_all_subj(:, i) = vs_allRuns;
    clear vs_allRuns   
end   

label = cell(1, length(dislike_all_subj(1).trial{1,1}(:,1)));

for k=1:length(dislike_all_subj(1).trial{1,1}(:,1))
    label{k} = num2str(k);
end

for k = 1:length(subjList) 
    dislike_all_subj(k).label = label;
end
    
size_dislike_all_subj = size(dislike_all_subj);
cfg = [];
dislike_all_subj_appended = ft_appenddata(cfg, dislike_all_subj(:,1), dislike_all_subj(:,2));

for k = 3:size_dislike_all_subj (2)
    dislike_all_subj_appended = ft_appenddata(cfg, dislike_all_subj_appended, dislike_all_subj(:,k)) ;   
end
clear dislike_all_subj

    
 %% like:

like_all_subj = struct('trial', [], 'time', []);

for i =  1:length(subjList) 
    load ([path2disc subjList(i).name '\runs_appended\virtsens\' cfg_virtsens '_', 'like_allRuns_' freq '.mat'], 'vs_allRuns');
    like_all_subj(:, i) = vs_allRuns;
    clear vs_allRuns   
end  
for k = 1:length(subjList) 
    like_all_subj(k).label = label;
end

size_like_all_subj = size(like_all_subj);
cfg = [];
like_all_subj_appended = ft_appenddata(cfg, like_all_subj(:,1), like_all_subj(:,2));

for k = 3:size_dislike_all_subj (2)
    like_all_subj_appended = ft_appenddata(cfg, like_all_subj_appended, like_all_subj(:,k)) ;   
end
clear like_all_subj

    
end

