function for_all_subj (subjpath, result_path, filename, freq, path2fig, algo, begin_time, end_time, perc_acc, Subject_num)

%% produziert Abbildungen für time series crossvalidation eines jeden Probanden 
   
    list = dir(subjpath);
    list(1:2)=[];
    like_vs_dislike_all = [];
    like_vs_dontcare_all = [];
    dislike_vs_dontcare_all = [];
    accuracy_counts = [];
    multibar = [];
 

 %% like vs. dislike
 
 
for m = 1: length(algo)
    for i = 1:length(list)   
           path2SVMresult = ([subjpath list(i).name result_path]); 
           [result] = kh_figure(path2SVMresult, filename, freq);
           like_vs_dislike_all.(algo{m}).(list(i).name).accuracy = squeeze(result.(lower(algo{m}))(1,2,:));
           like_vs_dislike_all.(algo{m}).(list(i).name).name = list(i).name;
    end

     
    for k = 1:numel(fieldnames(like_vs_dislike_all.(algo{m})))
        Accuracy_all.(algo{m})(k,:) = like_vs_dislike_all.(algo{m}).(list(k).name).accuracy;
        Name{k} = like_vs_dislike_all.(algo{m}).(list(k).name).name;
    end

    % Proband Nr. 5 hat kein like:
    Accuracy_all.(algo{m})(5,:)=[];

    for tt = 1:length(begin_time) 
       [accuracy_counts, multibar] = bar_time_interval(multibar, accuracy_counts, Accuracy_all, begin_time(1,tt),end_time(1, tt), algo, m, result, perc_acc, freq, path2fig );  
    end
    
    
end

fieldname = fields(multibar);
ylabels = ({'SVM', 'WeiRD', 'GNB'});
for op = 1:length(fields(multibar))
    figure
    bar(([sum(multibar.(fieldname{op}).svm); sum(multibar.(fieldname{op}).WeiRD); sum(multibar.(fieldname{op}).GNB)]));
    xticklabels(ylabels); 
    legend(multibar.(fieldname{op}).ylabel)
    legend boxoff
    grid on;
    ylabel(['n subjects']);
    title(['like_vs_dislike_' fieldname{op}]);
    savefig([path2fig 'like_vs_dislike_' fieldname{op} '_' freq '.fig'])
    close
end

end

function [accuracy_counts, multibar] = bar_time_interval(multibar, accuracy_counts, Accuracy_all, begin_time, end_time, algo, m, result, perc_acc, freq, path2fig)

   
    beg_sample = find(round(result.time{1,1}, 2) == (begin_time/1000)); 
    end_sample = find(round(result.time{1,1}, 2) == (end_time/1000));
    Accuracy_all_100 = 100*Accuracy_all.(algo{m})+50; 
    Accuracy_all_100 = Accuracy_all_100(:, beg_sample(1):end_sample(1));
    ind = zeros(size(Accuracy_all_100,1)+1, size(Accuracy_all_100,2));
    ind(end,:) = result.time{1,1}(beg_sample(1):end_sample(1));

    % adi_hist(Accuracy_all_100)

    % for i = 1:size(Accuracy_all_100,1)-1
    %     [~, y ] = find(Accuracy_all_100(i,:) >= 80);
    % %     z_ind = y(y>=beg_sample & y<=end_sample);
    %     ind(i, y) = 1; 
    %     clear y
    % end
    
    intervall = ['poststim_' num2str(begin_time) '_to_' num2str(end_time) '_ms'];

    for p=1:length(perc_acc)
        for i = 1:size(Accuracy_all_100,1)

            z_ind = Accuracy_all_100(i,:);
            z_ind = z_ind(z_ind >= perc_acc(1,p));
            if ~isempty(z_ind)
                accuracy_counts.(algo{m}).(intervall).(['perc_acc_from_' num2str(perc_acc(1,p))])(i) = true;
            else
                accuracy_counts.(algo{m}).(intervall).(['perc_acc_from_' num2str(perc_acc(1,p))])(i) = false;
            end
            clear z_ind
        end
    end


    fields = fieldnames(accuracy_counts.(algo{m}).(intervall));
    histo_count = accuracy_counts.(algo{m}).(intervall).(['perc_acc_from_' num2str(perc_acc(1,1))])';
    for op = 2: length(fieldnames(accuracy_counts.(algo{m}).(intervall)))
        if isfield(accuracy_counts.(algo{m}).(intervall), (fields{op}))
            histo_count = cat(2, histo_count, accuracy_counts.(algo{m}).(intervall).(fields{op})');
        end
    end

    for op = 1:length(perc_acc)
        leg{op} = ['ab ' num2str(perc_acc(1,op)) ' %'] ;
    end
    
    multibar.(intervall).(algo{m}) = histo_count;
    multibar.(intervall).ylabel = leg;
    
%     leg = categorical(leg);% ({'apples','pears','oranges'});
%     figure
%     bar(leg, sum(histo_count));
%     grid on;
%     xlabel(['accuracy [%]']);
%     ylabel(['n subjects']);
%     title(['like_vs_dislike_' (algo{m}) '_' intervall]);
%     
%     savefig ([path2fig 'like_vs_dislike_' (algo{m}) '_' intervall '_' freq '.fig']);
%     close
    % figure
    % subplot(2,2,1)       % add first plot in 2 x 2 grid
    % bar(sum(histo_count));        % line plot
    
    
%     clear leg
end



function [result] = kh_figure (path2fig, filename, freq)
   
    if exist ([path2fig filename freq '.mat'], 'file')
         load ([path2fig filename freq '.mat']);
    else
        result = [];
    end

end