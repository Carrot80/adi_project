function  [] = grandavg_sensorspace(subjectpath, dir_input, path2save, comp)

subject_values.like = [];
subject_values.dislike = [];

    
for ii = 1:length(subjectpath)
    
    load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep dir_input 'grandavg_like.mat']);
    load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep dir_input 'grandavg_dislike.mat']);
    load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep dir_input 'grandavg_dontcare.mat']);
    
    %% like
    for kk = 1:length(comp)
        time_beg = nearest(avg_like.time, comp(kk,1));
        time_end = nearest(avg_like.time, comp(kk,2));
        values_like(kk,:) = avg_like.avg(time_beg:time_end);
        
    end
    
    counter = 1;
 


end  