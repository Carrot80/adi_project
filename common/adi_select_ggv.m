  
function [vs_allRuns_trigger] = regroup_trials(vs_allRuns, balldesign) 

 ind = find(vs_allRuns.triggerlabel == balldesign.triggerchannel);
 vs_allRuns_trigger.trial = vs_allRuns.trial(ind);
 vs_allRuns_trigger.time = vs_allRuns.time(ind);
 vs_allRuns_trigger.balldesign = vs_allRuns.balldesign(ind);
 vs_allRuns_trigger.triggerlabel = vs_allRuns.triggerlabel(ind);
 vs_allRuns_trigger.response_label = vs_allRuns.response_label(ind);
 vs_allRuns_trigger.fsample = vs_allRuns.fsample;
 vs_allRuns_trigger.label = vs_allRuns.label;
 clearvars ind
 
 if 1 == isequal(balldesign.triggerchannel, 108)
     switch subject
         case 'nl_adi_04'
             vs_allRuns_trigger.response_label(1:length(vs_allRuns_trigger.trial)) = {'dontcare'};
     end
 elseif 1 == isequal(balldesign.triggerchannel, 4198)
     switch subject
         case 'nl_adi_17' % 
             ind = find(strcmp(vs_allRuns_trigger.response_label, 'like' ));
             vs_allRuns_trigger.trial(ind) = [];
             vs_allRuns_trigger.time(ind) = [];
             vs_allRuns_trigger.balldesign(ind) = [];
             vs_allRuns_trigger.triggerlabel(ind) = [];
             vs_allRuns_trigger.response_label(ind) = [];
         case  'nl_adi_16' , case 'nl_adi_06' 
             ind = find(strcmp(vs_allRuns_trigger.response_label, 'dislike' ));
             vs_allRuns_trigger.trial(ind) = [];
             vs_allRuns_trigger.time(ind) = [];
             vs_allRuns_trigger.balldesign(ind) = [];
             vs_allRuns_trigger.triggerlabel(ind) = [];
             vs_allRuns_trigger.response_label(ind) = [];
             
             
     end
     
     
     
 end
     
end


