function [delete_run] = boolean_operations(path2excelsheet, subject_list, balldesign)

response_tbl = readtable('U:\My Documents\MATLAB\eigene_Skripte\adi_project\excelsheets\subjects.xlsx');
% subj_responses.table = table2struct(response_tbl, 'ToScalar',true);

cleaned_table =  response_tbl;

ind_balldesign = find(strcmp(response_tbl.balldesign, balldesign));

figure
boxplot(response_tbl.like_allRuns(ind_balldesign))


figure
boxplot(response_tbl.like_allRuns, response_tbl.balldesign)

% boolean operation rules: Proband wird behalten, wenn er mind. 14 mal like
% angegeben hat, aber nur einmal dislike
% oder:
% wenn er mind. 14 mal dislike
% angegeben hat, aber nur einmal like


responses_allruns = [response_tbl.like_allRuns response_tbl.dislike_allRuns response_tbl.dontcare];

for k=1:length(responses_allruns)
    if responses_allruns(k,1) >= 14 && responses_allruns(k,2)<=1
        keep_subj_balldesign(k, 1) = 1;
        delete_sub.(response_tbl.Subject{k}).(response_tbl.balldesign{k}) = 0;
    elseif responses_allruns(k,2) >= 14 && responses_allruns(k,1)<=1
        keep_subj_balldesign(k, 1) = 1;
        delete_sub.(response_tbl.Subject{k}).(response_tbl.balldesign{k}) = 0;
    elseif responses_allruns(k,1) <=1
        delete_sub.(response_tbl.Subject{k}).(response_tbl.balldesign{k}) = 0;
        keep_subj_balldesign(k, 1) = 1;
    else
        keep_subj_balldesign(k, 1) = 0;  
        delete_balldesign.(response_tbl.Subject{k}).(response_tbl.balldesign{k})=1;
    end
end

% for k = 1:length(keep_subj_balldesign)   
%     if 0==keep_subj_balldesign(k)
%         delete_balldesign.(response_tbl.Subject{k}).(response_tbl.balldesign{k})=1;
% %     else
% %         delete_sub.(response_tbl.Subject{k}).(response_tbl.balldesign{k})=0;
%     end 
% end

responses_run1 = [response_tbl.like_run1 response_tbl.dislike_run1 response_tbl.dontcare_run1];
 
for k=1:length(response_tbl.like_run1)
    if response_tbl.like_run1(k) >= 6 %&& responses_run1(k,2)<2
        keep_run1(k, 1) = 1;
    elseif response_tbl.dislike_run1(k) >= 2 && response_tbl.like_run1(k) <=1
        keep_run1(k, 1) = 1;
    else
        keep_run1(k, 1) = 0;
    end

end

responses_run2 = [response_tbl.like_run2 response_tbl.dislike_run2 response_tbl.dontcare_run2];
 
for k=1:length(responses_run2)
    if response_tbl.like_run2(k) >= 6 
        keep_run2(k, 1) = 1;
    elseif response_tbl.dislike_run2(k) >= 2 &&  response_tbl.like_run2(k) <=1
        keep_run2(k, 1) = 1;
    else
        keep_run2(k, 1) = 0;
    end

end

responses_run3 = [response_tbl.like_run3 response_tbl.dislike_run3 response_tbl.dontcare_run3];
 
for k=1:length(responses_run3)
    if response_tbl.like_run3(k) >= 6 
        keep_run3(k, 1) = 1;
    elseif response_tbl.dislike_run3(k) >= 2 response_tbl.like_run3(k)<=1
        keep_run3(k, 1) = 1;
    else
        keep_run3(k, 1) = 0;
    end

end

for k = 1:length(keep_run1)   
    if 1==keep_subj_balldesign(k) && 0==keep_run1(k) 
        delete_run.(response_tbl.Subject{k}).run1.(response_tbl.balldesign{k})=1;
    elseif 1==keep_subj_balldesign(k) && 1==keep_run1(k)  
        delete_run.(response_tbl.Subject{k}).run1.(response_tbl.balldesign{k})=0;
    elseif 0==keep_subj_balldesign(k) && 1==keep_run1(k)
        delete_run.(response_tbl.Subject{k}).run1.(response_tbl.balldesign{k})=1;
    elseif 0==keep_subj_balldesign(k) && 0==keep_run1(k)
        delete_run.(response_tbl.Subject{k}).run1.(response_tbl.balldesign{k})=1;
    end 
    
     if 1==keep_subj_balldesign(k) && 0==keep_run2(k)  
        delete_run.(response_tbl.Subject{k}).run2.(response_tbl.balldesign{k})=1;
     elseif 1==keep_subj_balldesign(k) && 1==keep_run2(k)  
        delete_run.(response_tbl.Subject{k}).run2.(response_tbl.balldesign{k})=0;
     elseif 0==keep_subj_balldesign(k) && 1==keep_run2(k)  
        delete_run.(response_tbl.Subject{k}).run2.(response_tbl.balldesign{k})=1;
     elseif 0==keep_subj_balldesign(k) && 0==keep_run2(k)
        delete_run.(response_tbl.Subject{k}).run2.(response_tbl.balldesign{k})=1;
     end    
    
     if 1==keep_subj_balldesign(k) && 0==keep_run3(k)  
        delete_run.(response_tbl.Subject{k}).run3.(response_tbl.balldesign{k})=1;
     elseif 1==keep_subj_balldesign(k) && 1==keep_run3(k)  
        delete_run.(response_tbl.Subject{k}).run3.(response_tbl.balldesign{k})=0;
     elseif 0==keep_subj_balldesign(k) && 1==keep_run3(k)  
        delete_run.(response_tbl.Subject{k}).run3.(response_tbl.balldesign{k})=1;
     elseif 0==keep_subj_balldesign(k) && 0==keep_run3(k)
        delete_run.(response_tbl.Subject{k}).run3.(response_tbl.balldesign{k})=1;  
        
     end      
end

delete_run.readme = 'if delete_run == 1,then run will be deleted; if delete_run == 0, then run will be kept';




end