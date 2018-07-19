function for_all_subj (subjpath, condition, freq, time, path2fig)

%% produziert Abbildungen für time series crossvalidation eines jeden Probanden 
   
    list = dir(subjpath);
    list(1:2)=[];
    tlike_vs_dislike_all = [];
%     tlike_vs_dontcare_all=[];
%     tdislike_vs_dontcare_all =[];

 

 %% like vs. dislike

    for i = 1:length(list)   
           path2SVMresult = ([subjpath list(i).name '\MEG_analysis\noisereduced\1_95Hz\04_timelockstatistics_interp\']); 
           [tlike_vs_dislike] = kh_figure(path2SVMresult, condition, freq);
           tlike_vs_dislike_all.(list(i).name).accuracy = tlike_vs_dislike.Condition1vs2.Accuracy;
           tlike_vs_dislike_all.(list(i).name).binomial = tlike_vs_dislike.Condition1vs2.Binominal;
           tlike_vs_dislike_all.(list(i).name).binomial_sig = find (tlike_vs_dislike.Condition1vs2.Binominal<=0.05);   
           tlike_vs_dislike_all.(list(i).name).name = list(i).name;
    end
    
    Binomial_sig_all = nan(numel(fieldnames(tlike_vs_dislike_all)), length(tlike_vs_dislike_all.(list(i).name).accuracy)) ;
    
    for k = 1:numel(fieldnames(tlike_vs_dislike_all))
        if isempty (tlike_vs_dislike_all.(list(k).name).binomial_sig)
            Binomial_sig_all(k,:)= NaN;
        else
            Binomial_sig_all(k,tlike_vs_dislike_all.(list(k).name).binomial_sig ) = tlike_vs_dislike_all.(list(k).name).binomial(tlike_vs_dislike_all.(list(k).name).binomial_sig); 
        end
        Accuracy_all(k,:) = tlike_vs_dislike_all.(list(k).name).accuracy;
        Name{k} = tlike_vs_dislike_all.(list(k).name).name;
    end
    
 figure   
 plot (time, Accuracy_all)
%  hold on
%  plot(time, Binomial_sig_all, '+')     
%  legend(Name)
 title(['tlike_vs_tdislike ' freq]);
 ylim([0 1])
 xlabel('time');
 ylabel ('accuracy');
 set(findall(gca, 'Type', 'Line'),'LineWidth', 1);
 savefig([path2fig 'tlike_vs_dislike_' freq '_old_fig.fig'])
 close 
     
     %% plot figure with mean and std with shades:
     
    mean_Accuracy = mean(Accuracy_all);
    mean(Accuracy_all(:,1))
    std_Accuracy = std(Accuracy_all);
    std_Accuracy_pos = mean_Accuracy + std_Accuracy;
    std_Accuracy_neg = mean_Accuracy - std_Accuracy;  
    max_Accuracy = max(Accuracy_all); 
    min_Accuracy = min(Accuracy_all);
    
    figure
    [ph, msg] = jbfill(time, std_Accuracy_pos, std_Accuracy_neg, 'green', [0.1, 0.1, 0.1],1, 0.2)
    hold on
    plot(time, mean_Accuracy, 'k', 'LineWidth', 1) 
    hold on
    [ph, msg] = jbfill(time, max_Accuracy, std_Accuracy_pos, [0.1, 0.5, 0.5], [0.1, 0.1, 0.1],1, 0.2)
    hold on
    [ph, msg] = jbfill(time, std_Accuracy_neg, min_Accuracy, [0.1, 0.5, 0.5], [0.1, 0.1, 0.1],1, 0.2)
    grid on
    ylim([0.2 1])
    xlabel('time');
    ylabel ('Accuracy of Classification');
    title(['like vs dislike' '_' freq])
    savefig([path2fig 'tlike_vs_dislike_' freq '.fig'])
    close  
 %     [ph,msg]=jbfill(time,std_Accuracy_pos,std_Accuracy_neg,rand(1,3),rand(1,3),0,rand(1,1))
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,0.5)    
         



%      %% like vs dontcare
% figure;
%   for i=1:length(list)  
%            
%            [tlike_vs_dontcare] = kh_figure(list(i).name, filter, 'tlike_vs_tdontcare.mat');
%             if ~isempty(tlike_vs_dontcare)
%                tlike_vs_dontcare_all.(list(i).name).accuracy=tlike_vs_dontcare.tlike_vs_tdontcare.MEG.Accuracy;
%                tlike_vs_dontcare_all.(list(i).name).binomial=tlike_vs_dontcare.tlike_vs_tdontcare.MEG.Binominal;
%                tlike_vs_dontcare_all.(list(i).name).binomial_sig = find (tlike_vs_dontcare.tlike_vs_tdontcare.MEG.Binominal<=0.05);
%                tlike_vs_dontcare_all.(list(i).name).name = list(i).name;
%             end
%   end         
%             
%     Binomial_sig_all = nan( numel(fieldnames(tlike_vs_dontcare_all)), length(tlike_vs_dontcare_all.(list(i).name).accuracy)) ;
%     
%     list_dontcare = fieldnames(tlike_vs_dontcare_all); 
% 
%     
%     for k=1:numel(fieldnames(tlike_vs_dontcare_all))
%             if isempty (tlike_vs_dontcare_all.(list_dontcare{k}).binomial_sig)
%             Binomial_sig_all(k,:)= NaN;
%         else 
%             Binomial_sig_all(k,tlike_vs_dontcare_all.(list_dontcare{k}).binomial_sig ) = tlike_vs_dontcare_all.(list_dontcare{k}).binomial(tlike_vs_dontcare_all.(list_dontcare{k}).binomial_sig); 
%         end
%         Accuracy_all(k,:) = tlike_vs_dontcare_all.(list_dontcare{k}).accuracy;
%         Name_dontcare{k}=tlike_vs_dontcare_all.(list_dontcare{k}).name;
%     end  
%             
%      graph_1 = plot (time, Accuracy_all);
%      hold on
%      plot(time, Binomial_sig_all, '+')     ;
%      legend(Name_dontcare);
%      title('like_vs_dontcare');
%      ylim([0 1])
%      xlabel('time');
%      ylabel ('accuracy, p-value');
%      set(findall(gca, 'Type', 'Line'),'LineWidth', 1.5);
% %      set(graph_1, 'LineWidth', 1.5);
%      
%      clear Accuracy_all Binomial_sig_all  tlike_vs_dontcare_all
% 
%    %% dislike vs. dontcare
%    
%   figure
%   for i=1:length(list)   
%            [tdislike_vs_dontcare] = kh_figure(list(i).name, filter, 'tdislike_vs_tdontcare.mat')
%             if ~isempty(tdislike_vs_dontcare)
%                tdislike_vs_dontcare_all.(list(i).name).accuracy=tdislike_vs_dontcare.tdislike_vs_tdontcare.MEG.Accuracy;
%                tdislike_vs_dontcare_all.(list(i).name).binomial=tdislike_vs_dontcare.tdislike_vs_tdontcare.MEG.Binominal;
%                tdislike_vs_dontcare_all.(list(i).name).binomial_sig = find (tdislike_vs_dontcare.tdislike_vs_tdontcare.MEG.Binominal<=0.05);
%                tdislike_vs_dontcare_all.(list(i).name).name = list(i).name;
%             end
%   end
%   
%     Binomial_sig_all = nan( numel(fieldnames(tdislike_vs_dontcare_all)), length(tdislike_vs_dontcare_all.(list(i).name).accuracy)) ;
%     list_dontcare = fieldnames(tdislike_vs_dontcare_all); 
%   
%     for k=1:numel(fieldnames(tdislike_vs_dontcare_all))
%         if isempty (tdislike_vs_dontcare_all.(list_dontcare{k}).binomial_sig)
%             Binomial_sig_all(k,:)= NaN;
%         else
%             Binomial_sig_all(k,tdislike_vs_dontcare_all.(list_dontcare{k}).binomial_sig ) = tdislike_vs_dontcare_all.(list_dontcare{k}).binomial(tdislike_vs_dontcare_all.(list_dontcare{k}).binomial_sig); 
%         end
%         Accuracy_all(k,:) = tdislike_vs_dontcare_all.(list_dontcare{k}).accuracy;
%         Name_dontcare{k}=tdislike_vs_dontcare_all.(list_dontcare{k}).name;
%     end  
%             
%      graph_2 = plot (time, Accuracy_all);
%      hold on
%      plot(time, Binomial_sig_all, '+');     
%      legend(Name_dontcare);
%      title('dislike_vs_dontcare');
%      ylim([0 1]);
%      xlabel('time');
%      ylabel ('accuracy, p-value');
%      set(findall(gca, 'Type', 'Line'),'LineWidth', 1.5);
% %      set(graph_1, 'LineWidth', 1.5);
%  
end


function [condition1_vs_condition2] = kh_figure (path2fig, condition, freq)
   
    if exist ([path2fig filesep condition '_' freq '.mat'], 'file')
         condition1_vs_condition2 = load ([path2fig filesep condition '_' freq '.mat']);
    else
        condition1_vs_condition2 = [];
    end

end