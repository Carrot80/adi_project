
function adi_append_group_MEG (fieldtripPath, outPath, bpfreq)

subjList = dir( fieldtripPath );
subjList(1:2)=[]; % remove the first two ./..
    
%%   dislike:
outFile_dislike = strcat(outPath, 'dislike_group_', bpfreq, '.mat');
if ~exist (outFile_dislike, 'file')
    dislike_all_subj = struct('label', [] , 'dimord', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', []);
 
    for i =  1:length(subjList) 
        [dislike_all_subj]  =  kh_appendalldata( strcat(fieldtripPath, subjList(i,1).name, '\MEG_analysis\noisereduced\', '1_95Hz', '\03_appended_data\'), subjList(i,1).name, bpfreq, dislike_all_subj, 'dislike_allRuns', i)
    end   
    
    save(outFile_dislike , 'dislike_all_subj', '-v7.3');
       
%     size_dislike_all_subj=size(dislike_all_subj);
%     cfg=[];
%     dislike_all_subj_appended = ft_appenddata(cfg, dislike_all_subj(:,1), dislike_all_subj(:,2));
%     
%     for k = 3:size_dislike_all_subj (2)
%         dislike_all_subj_appended = ft_appenddata(cfg, dislike_all_subj_appended, dislike_all_subj(:,k)) ;   
%     end
%     save(outFile_dislike , 'dislike_all_subj_appended', '-v7.3');
%     clear dislike_all_subj_appended
end
    
 %% like:
outFile_like = strcat(outPath, 'like_group_', bpfreq, '.mat');
if ~exist (outFile_like, 'file')
    like_all_subj = struct('label', [] , 'dimord', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', []);

    for i =  1:length(subjList) 
        [like_all_subj]  =  kh_appendalldata( strcat(fieldtripPath, subjList(i,1).name, '\MEG_analysis\noisereduced\', '1_95Hz', '\03_appended_data\'), subjList(i,1).name, bpfreq, like_all_subj, 'like_allRuns', i);
    end   
    save ([outPath, 'like_group_', bpfreq], 'like_all_subj', '-v7.3')
    
%     size_like_all_subj=size(like_all_subj);
%     cfg=[];
%     like_all_subj_appended = ft_appenddata(cfg, like_all_subj(:,1), like_all_subj(:,2));
%     
%     for k = 3:size_like_all_subj (2)
%         like_all_subj_appended = ft_appenddata(cfg, like_all_subj_appended, like_all_subj(:,k)) ;   
%     end
%     save(outFile_like , 'like_all_subj_appended', '-v7.3');
%     clear like_all_subj_appended
end
    
 
    
end


function [condition_all_subj]=kh_appendalldata ( subjFolder, subj, bpfreq, condition_all_subj, name_condition, i)
 
condition  =  load(strcat(subjFolder, strcat(name_condition, '_', bpfreq, '.mat')));
if ~isfield(condition.(name_condition), 'sampleinfo')
    condition.(name_condition) = setfield(condition.(name_condition), 'sampleinfo', []);
end
if ~isfield(condition.(name_condition), 'dimord')
    condition.(name_condition) = setfield(condition.(name_condition), 'dimord', 'rpt_chan_time');
end
condition.(name_condition) = orderfields (condition.(name_condition), {'label' , 'dimord', 'sampleinfo', 'trial', 'time', 'cfg' });
condition_all_subj(i) = condition.(name_condition);
condition.(name_condition).subject = subj

end
