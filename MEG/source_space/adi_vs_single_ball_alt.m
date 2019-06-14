
function  [virtsens_all_subj] = adi_appendvirtsensMEG(virtsens_all_subj, path2virtsens, dir_out, subject, pattern, trigger, freq, balldesign, delete_runs, subjnum)

dir_runs = dir(path2virtsens);
dir_runs(1:2)=[];

for i = 1:length(delete_runs)
    if 1==strcmp(delete_runs{i,1}, subject)
        dir_runs(delete_runs{i,2},:)=[];
    end
end

for k = 1:length(dir_runs)
    virtsens(k)= load ([path2virtsens  dir_runs(k).name filesep 'virtsens_' freq '.mat'], 'virtsens');
end

switch length(dir_runs)
    case 3
        vs_allRuns.trial = [virtsens(1).virtsens.trial  virtsens(2).virtsens.trial virtsens(3).virtsens.trial];
        vs_allRuns.time = [virtsens(1).virtsens.time  virtsens(2).virtsens.time virtsens(3).virtsens.time];
        vs_allRuns.response = [virtsens(1).virtsens.response virtsens(2).virtsens.response virtsens(3).virtsens.response];
        vs_allRuns.response_label = [virtsens(1).virtsens.response_label virtsens(2).virtsens.response_label virtsens(3).virtsens.response_label];
        vs_allRuns.balldesign = [virtsens(1).virtsens.balldesign virtsens(2).virtsens.balldesign virtsens(3).virtsens.balldesign];
        
    case 2
        vs_allRuns.trial = [virtsens(1).virtsens.trial  virtsens(2).virtsens.trial];
        vs_allRuns.time = [virtsens(1).virtsens.time  virtsens(2).virtsens.time];
        vs_allRuns.response = [virtsens(1).virtsens.response virtsens(2).virtsens.response];
        vs_allRuns.response_label = [virtsens(1).virtsens.response_label virtsens(2).virtsens.response_label];
        vs_allRuns.balldesign = [virtsens(1).virtsens.balldesign virtsens(2).virtsens.balldesign];
    case 1
        vs_allRuns.trial = virtsens(1).virtsens.trial;
        vs_allRuns.time = virtsens(1).virtsens.time;
        vs_allRuns.response = virtsens(1).virtsens.response;
        vs_allRuns.response_label = virtsens(1).virtsens.response_label;
        vs_allRuns.balldesign = virtsens(1).virtsens.balldesign;
end
    vs_allRuns.label = virtsens(1).virtsens.label;
    vs_allRuns.fsample = virtsens(1).virtsens.fsample;
    
    clearvars virtsens
    
    for k=1:length(trigger.labels_short)
       temp(k)=strcmp(trigger.labels_short{k}, balldesign);
    end
    ind=find(temp);
    balldesign_long = trigger.labels{ind};
    clear ind
    for k=1:length(vs_allRuns.trial)
        ind(k) = strcmp(vs_allRuns.balldesign{k}, balldesign_long);
    end
  
vs_allRuns_ball.trial =  vs_allRuns.trial(find(ind));
vs_allRuns_ball.time =  vs_allRuns.time(find(ind));
vs_allRuns_ball.response =  vs_allRuns.response(find(ind));
vs_allRuns_ball.response_label =  vs_allRuns.response_label(find(ind));
vs_allRuns_ball.balldesign =  vs_allRuns.balldesign(find(ind));
vs_allRuns_ball.subject = subject;

virtsens_all_subj(subjnum)=vs_allRuns_ball;

end

     
