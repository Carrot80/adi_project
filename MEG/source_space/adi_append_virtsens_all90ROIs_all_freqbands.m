
 function [session]= main(sessions, path2data, outPath_extdisc, freqname, conditions, meanmax)

 [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, conditions, num2str(1), meanmax);
 [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, conditions, num2str(2), meanmax);
 [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, conditions, num2str(3), meanmax);
 [session] = kh_concatenate_runs(sessions);

 
 end
 
 
 function [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, conditions, run, meanmax)
  
 for m = 1:size(freqname,1)
  
  fileName = (['vs_' meanmax '*' freqname{m} '.mat']);
  files = dir(fullfile([path2data '90ROIs_plus_freqbands\run' run filesep fileName]));
  size_files = size(files);
    
  for i = 1:(size_files(1,1))
      virtsens = load ([path2data '90ROIs_plus_freqbands\run' run filesep files(i).name]);
        if 1 == strcmp(freqname{m}, 'bp1-45Hz')
            freqname{m} = 'bp1_45Hz';
        end
        if 1 == contains(files(i).name, '_dislike_')
            condition = 2;
            fn = [meanmax '_all_rois_dislike' ];
            for k = 1:size(virtsens.(fn).trial,2)
                for p = 1:90 %size(virtsens.(fn).trial(:,k),1)
                    data.dislike.(freqname{m})(k,p,:) = virtsens.(fn).trial{p,k};
                end 
            end
            labels.dislike.(freqname{m}) = 2*ones(1,size(virtsens.(fn).trial,2));
        elseif 1 == contains(files(i).name, '_like_')
            condition = 1;
            fn = [meanmax '_all_rois_like' ];
            for k = 1:size(virtsens.(fn).trial,2)
                for p = 1:90%size(virtsens.(fn).trial(:,k),1)
                    data.like.(freqname{m})(k,p,:) = virtsens.(fn).trial{p,k};
                end 
            end
            labels.like.(freqname{m}) = ones(1,size(virtsens.(fn).trial,2));
             
        elseif 1 == contains(files(i).name, '_dontcare_')
            condition = 3;
             fn = [meanmax '_all_rois_dontcare' ];
             for k = 1:size(virtsens.(fn).trial,2)
                 for p = 1:90%size(virtsens.(fn).trial(:,k),1)
                    data.dontcare.(freqname{m})(k,p,:) = virtsens.(fn).trial{p,k};
                 end 
             end
             labels.dontcare.(freqname{m}) = 3*ones(1,size(virtsens.(fn).trial,2));
        end 
          
%         size_sessions_data = size(sessions(str2double(run)).data, 1);
%         trl_num = size(virtsens.(fn).trial,2);
%         start = 1+size_sessions_data;
%         end_ = start + trl_num -1;
%         for n=start:end_
%             sessions(str2double(run)).data(n,:,:) = virtsens.(fn).trial{1,n-start+1};
%         end

  end
  
num_conditions = size(fields(data),1);
fieldnames =   fields(data);  
switch num_conditions
    case 2
        sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).data = data.(fieldnames{1}).(freqname{m});
        sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).labels = labels.(fieldnames{1}).(freqname{m});
        sessions.(fieldnames{2}).(freqname{m}).(['run_' run]).data = data.(fieldnames{2}).(freqname{m});
        sessions.(fieldnames{2}).(freqname{m}).(['run_' run]).labels = labels.(fieldnames{2}).(freqname{m});
    case 3
        sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).data = data.(fieldnames{1}).(freqname{m});
        sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).labels = labels.(fieldnames{1}).(freqname{m});
        sessions.(fieldnames{2}).(freqname{m}).(['run_' run]).data = data.(fieldnames{2}).(freqname{m});
        sessions.(fieldnames{2}).(freqname{m}).(['run_' run]).labels = labels.(fieldnames{2}).(freqname{m});
        sessions.(fieldnames{3}).(freqname{m}).(['run_' run]).data = data.(fieldnames{3}).(freqname{m});
        sessions.(fieldnames{3}).(freqname{m}).(['run_' run]).labels = labels.(fieldnames{3}).(freqname{m});

end

sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).tissuelabel = virtsens.(fn).tissuelabel; 
sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).time = virtsens.(fn).time; 
sessions.(fieldnames{1}).(freqname{m}).(['run_' run]).meanmax = meanmax;
%   sessions(str2double(run)).labels(start:end_) = condition.* ones(1, length(virtsens.trial));
  
 end   
end

function  [session] = kh_concatenate_runs(sessions)

num_cond = size(fields(sessions),1);
fieldnames_conditions = fields(sessions);  
num_freqbands = size(fields(sessions.(fieldnames_conditions{1})),1);
fieldnames_freqbands = fields(sessions.(fieldnames_conditions{1}));
num_runs = size(fields(sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1})),1);
fieldnames_runs = fields(sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1}));

switch num_cond
    case 3 
        for i = 1:num_freqbands
            data_like.(fieldnames_freqbands{i}).data = sessions.like.(fieldnames_freqbands{i,1}).(fieldnames_runs{1}).data;
            data_dislike.(fieldnames_freqbands{i}).data = sessions.dislike.(fieldnames_freqbands{i,1}).(fieldnames_runs{1}).data;
            data_dontcare.(fieldnames_freqbands{i}).data = sessions.dontcare.(fieldnames_freqbands{i,1}).(fieldnames_runs{1}).data;
            labels_like = sessions.like.(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels;
            labels_dislike = sessions.dislike.(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels;
            labels_dontcare = sessions.dontcare.(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels;

            for k = 2:num_runs
                if isfield (sessions.like.(fieldnames_freqbands{i,1}), (fieldnames_runs{1}))
                    data_like.(fieldnames_freqbands{i}).data = cat(1, data_like.(fieldnames_freqbands{i}).data, sessions.like.(fieldnames_freqbands{i}).(fieldnames_runs{k}).data);
                    data_dislike.(fieldnames_freqbands{i}).data = cat(1, data_dislike.(fieldnames_freqbands{i}).data, sessions.dislike.(fieldnames_freqbands{i}).(fieldnames_runs{k}).data);
                    data_dontcare.(fieldnames_freqbands{i}).data = cat(1, data_dontcare.(fieldnames_freqbands{i}).data, sessions.dontcare.(fieldnames_freqbands{i}).(fieldnames_runs{k}).data);
                    labels_like = cat(2, labels_like, sessions.like.(fieldnames_freqbands{1}).(fieldnames_runs{k}).labels);
                    labels_dislike = cat(2, labels_dislike, sessions.dislike.(fieldnames_freqbands{1}).(fieldnames_runs{k}).labels); 
                    labels_dontcare = cat(2, labels_dontcare, sessions.dontcare.(fieldnames_freqbands{1}).(fieldnames_runs{k}).labels); 
                end
            end
        end
        labels_like = cat(2, sessions.like.(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels, sessions.like.(fieldnames_freqbands{1}).(fieldnames_runs{2}).labels, sessions.like.(fieldnames_freqbands{1}).(fieldnames_runs{3}).labels);
        labels_dislike = cat(2, sessions.dislike.(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels, sessions.dislike.(fieldnames_freqbands{1}).(fieldnames_runs{2}).labels, sessions.dislike.(fieldnames_freqbands{1}).(fieldnames_runs{3}).labels);
        labels_dontcare = cat(2, sessions.dontcare.(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels, sessions.dontcare.(fieldnames_freqbands{1}).(fieldnames_runs{2}).labels, sessions.dontcare.(fieldnames_freqbands{1}).(fieldnames_runs{3}).labels);
        data_like_all_freqbands = data_like.(fieldnames_freqbands{1}).data;
        data_dislike_all_freqbands = data_dislike.(fieldnames_freqbands{1}).data;
        data_dontcare_all_freqbands = data_dontcare.(fieldnames_freqbands{1}).data;

        for k = 2:num_freqbands
            if isfield(data_like,(fieldnames_freqbands{k}))
                data_like_all_freqbands = cat(2, data_like_all_freqbands, data_like.(fieldnames_freqbands{k}).data);
                data_dislike_all_freqbands = cat(2, data_dislike_all_freqbands, data_dislike.(fieldnames_freqbands{k}).data);
                data_dontcare_all_freqbands = cat(2, data_dontcare_all_freqbands, data_dontcare.(fieldnames_freqbands{k}).data);
            end
        end
        
        session.labels = cat(2, labels_like, labels_dislike, labels_dontcare); 
        session.data = cat(1, data_like_all_freqbands, data_dislike_all_freqbands, data_dontcare_all_freqbands);           
 %%       
   case 2 
       
         for i = 1:num_freqbands
            data_condition_1.(fieldnames_freqbands{i}).data = sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{i,1}).(fieldnames_runs{1}).data;
            data_condition_2.(fieldnames_freqbands{i}).data = sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{i,1}).(fieldnames_runs{1}).data;
            labels_condition_1 = sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels;
            labels_condition_2 = sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels;
           
            for k = 2:num_runs
                if isfield (sessions.like.(fieldnames_freqbands{i,1}), (fieldnames_runs{1}))
                    data_condition_1.(fieldnames_freqbands{i}).data = cat(1, data_condition_1.(fieldnames_freqbands{i}).data, sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{i}).(fieldnames_runs{k}).data);
                    data_condition_2.(fieldnames_freqbands{i}).data = cat(1, data_condition_2.(fieldnames_freqbands{i}).data, sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{i}).(fieldnames_runs{k}).data);
                    labels_condition_1 = cat(2, labels_condition_1, sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1}).(fieldnames_runs{k}).labels);
                    labels_condition_2 = cat(2, labels_condition_2, sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{1}).(fieldnames_runs{k}).labels); 
                end
            end
          end
        labels_condition_1 = cat(2, sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels, sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1}).(fieldnames_runs{2}).labels, sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1}).(fieldnames_runs{3}).labels);
        labels_condition_2 = cat(2, sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).labels, sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{1}).(fieldnames_runs{2}).labels, sessions.(fieldnames_conditions{2}).(fieldnames_freqbands{1}).(fieldnames_runs{3}).labels);
   
        data_condition_1_all_freqbands = data_condition_1.(fieldnames_freqbands{1}).data;
        data_condition_2_all_freqbands = data_condition_2.(fieldnames_freqbands{1}).data;
                 
        for k = 2:num_freqbands
            if isfield(data_condition_1,(fieldnames_freqbands{k}))
                data_condition_1_all_freqbands = cat(2, data_condition_1_all_freqbands, data_condition_1.(fieldnames_freqbands{k}).data);
                data_condition_2_all_freqbands = cat(2, data_condition_2_all_freqbands, data_condition_2.(fieldnames_freqbands{k}).data);
            end
        end
        session.labels = cat(2, labels_condition_1, labels_condition_2); 
        session.data = cat(1, data_condition_1_all_freqbands, data_condition_2_all_freqbands);
      
end
   
 session.tissuelabel = sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).tissuelabel;   
 session.time = sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).time;       
 session.meanmax = sessions.(fieldnames_conditions{1}).(fieldnames_freqbands{1,1}).(fieldnames_runs{1}).meanmax;      
 session.conditions = fieldnames_conditions;
 session.labels_readme = {'like = 1'; 'dislike = 2'; 'dontcare = 3'};
 session.ROIs = repmat(session.tissuelabel,1,num_freqbands)';
    num_ROIs = size(session.tissuelabel,2);
    begin = 1;
    k=1;
    for i = 1:num_freqbands
        session.ROIs(begin:num_ROIs*k,2) = repmat(fieldnames_freqbands(i),1, num_ROIs)';
        begin = begin+num_ROIs;
        k=k+1;
    end


end





     