function for_all_subj (subjpath, filename, freq, path2fig, algo)

%% produziert Abbildungen für time series crossvalidation eines jeden Probanden 
   
    list = dir(subjpath);
    list(1:2)=[];
    like_vs_dislike_all = [];
    like_vs_dontcare_all = [];
    dislike_vs_dontcare_all = [];

 

 %% like vs. dislike
for m = 1: length(algo)
    for i = 1:length(list)   
           path2SVMresult = ([subjpath list(i).name '\MEG\sourcespace\noRois\Guggenmos_decoding_results\']); 
           [result] = kh_figure(path2SVMresult, filename, freq);
           like_vs_dislike_all.(algo{m}).(list(i).name).accuracy = squeeze(result.(lower(algo{m}))(1,2,:));
           like_vs_dislike_all.(algo{m}).(list(i).name).name = list(i).name;
    end
    
%     Binomial_sig_all = nan(numel(fieldnames(tlike_vs_dislike_all)), length(tlike_vs_dislike_all.(list(i).name).accuracy)) ;
%     
    for k = 1:numel(fieldnames(like_vs_dislike_all.(algo{m})))
        Accuracy_all.(algo{m})(k,:) = like_vs_dislike_all.(algo{m}).(list(k).name).accuracy;
        Name{k} = like_vs_dislike_all.(algo{m}).(list(k).name).name;
    end
    
    % Proband nr. 5 (nl_adi_08) hat keine likes:
    Accuracy_all.(algo{m})(5,:) = [];
    figure; axis tight; hold on;   
    plot (result.time{1,1}, 100*Accuracy_all.(algo{m})+50)
    plot([-0.5 1], [50 50], 'k-')
%     set(findall(gca, 'Type', 'Line'),'LineWidth', 1);
    xlim([-0.5 1])
    xlabel('Time [s]');
    ylabel(['Classification accuracy [%]']);
    title([(algo{m}) ' 14 subjects'])
    savefig([path2fig 'decoding_result_like_vs_dislike_' freq '_14_subj_' (algo{m}) '.fig'])
    close
    

     
     %% plot figure with mean and std with shades:
     
    mean_Accuracy.(algo{m}) = mean(Accuracy_all.(algo{m}));
    mean(Accuracy_all.(algo{m})(:,1))
    std_Accuracy.(algo{m}) = std(Accuracy_all.(algo{m}));
    std_Accuracy_pos.(algo{m}) = mean_Accuracy.(algo{m}) + std_Accuracy.(algo{m});
    std_Accuracy_neg.(algo{m}) = mean_Accuracy.(algo{m}) - std_Accuracy.(algo{m});  
    max_Accuracy.(algo{m}) = max(Accuracy_all.(algo{m})); 
    min_Accuracy.(algo{m}) = min(Accuracy_all.(algo{m}));
    
    figure; axis tight;
    [ph, msg] = jbfill(result.time{1,1}, 100*std_Accuracy_pos.(algo{m})+50, 100*std_Accuracy_neg.(algo{m})+50, [0.1, 0.5, 0.5], [0.1, 0.1, 0.1],1, 0.2);
    hold on
    plot(result.time{1,1}, 100*mean_Accuracy.(algo{m})+50, 'k', 'LineWidth', 1) 
    hold on
    [ph, msg] = jbfill(result.time{1,1}, 100*max_Accuracy.(algo{m})+50, 100*std_Accuracy_pos.(algo{m})+50, 'magenta', [0.1, 0.1, 0.1],1, 0.2);
    hold on
    [ph, msg] = jbfill(result.time{1,1}, 100*std_Accuracy_neg.(algo{m})+50, 100*min_Accuracy.(algo{m})+50, 'magenta', [0.1, 0.1, 0.1],1, 0.2);
    grid on
    ylim([0 100])
    xlabel('time [s]');
    ylabel ('Accuracy of Classification [%]');
    title(['like vs dislike_' (algo{m}) '_' freq])
    savefig([path2fig (algo{m}) '_tlike_vs_dislike_' freq '.fig'])
    close  
 %     [ph,msg]=jbfill(time,std_Accuracy_pos,std_Accuracy_neg,rand(1,3),rand(1,3),0,rand(1,1))
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,0.5)    



end


end


function [result] = kh_figure (path2fig, filename, freq)
   
    if exist ([path2fig filename freq '.mat'], 'file')
         load ([path2fig filename freq '.mat']);
    else
        result = [];
    end

end