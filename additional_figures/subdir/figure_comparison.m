%%

function figure_comparison(tlike_vs_dislike_all, path, freq, MEEG)
subjects = fieldnames(tlike_vs_dislike_all);
data_freq = [];
fileName = (['like_vs_dislike_allRuns_', freq ]);
for i = 1:numel(fieldnames(tlike_vs_dislike_all))
    data_freq.accuracy(i,:) = tlike_vs_dislike_all.(subjects{i}).(fileName).accuracy
    data_freq.binomial(i,:) = tlike_vs_dislike_all.(subjects{i}).(fileName).binomial
    data_freq.binomial_sig(i,:) = tlike_vs_dislike_all.(subjects{i}).(fileName).binomial
end
figure
time = tlike_vs_dislike_all.(subjects{1}).(fileName).latency(1,:);
plot(time, data_freq.accuracy)
title(['tlike_vs_tdislike_', freq, '_', MEEG]);
ylim([0 1])
xlabel('time');
ylabel ('accuracy, p-value'); 
savefig([path, filesep, fileName, '.fig']);
fig = ([path, filesep, fileName]);
print('-dpng', fig); 
% 
% Binomial_sig_all = nan( numel(fieldnames(tlike_vs_dislike_all.(subj))), length(tlike_vs_dislike_all.(list(i).name).accuracy)) ;
%     
%     for k=1:numel(fieldnames(tlike_vs_dislike_all))
%         if isempty (tlike_vs_dislike_all.(list(k).name).binomial_sig)
%             Binomial_sig_all(k,:)= NaN;
%         else
%             Binomial_sig_all(k,tlike_vs_dislike_all.(list(k).name).binomial_sig ) = tlike_vs_dislike_all.(list(k).name).binomial(tlike_vs_dislike_all.(list(k).name).binomial_sig); 
%         end
%         Accuracy_all(k,:) = tlike_vs_dislike_all.(list(k).name).accuracy;
%         Name{k}=tlike_vs_dislike_all.(list(k).name).name;
%     end
%         
%         
%      plot (time, Accuracy_all)
%      hold on
%      plot(time, Binomial_sig_all, '+')     
%      legend(Name)
%      title('tlike_vs_tdislike');
%      ylim([0 1])
%      xlabel('time');
%      ylabel ('accuracy, p-value');
%      set(findall(gca, 'Type', 'Line'),'LineWidth', 1.5);
%      clear Accuracy_all Binomial_sig_all tlike_vs_dislike_all Name
% 
%      %% like vs dontcare
% 
%   for k=1:length(list)  
%            
%            [tlike_vs_dontcare] = kh_figure(list(k).name, filter, 'tlike_vs_tdontcare.mat');
%             if ~isempty(tlike_vs_dontcare)
%                tlike_vs_dontcare_all.(list(k).name).accuracy=tlike_vs_dontcare.tlike_vs_tdontcare.MEG.Accuracy;
%                tlike_vs_dontcare_all.(list(k).name).binomial=tlike_vs_dontcare.tlike_vs_tdontcare.MEG.Binominal;
%                tlike_vs_dontcare_all.(list(k).name).binomial_sig = find (tlike_vs_dontcare.tlike_vs_tdontcare.MEG.Binominal<=0.05);
%                tlike_vs_dontcare_all.(list(k).name).name = list(k).name;
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
%   for k=1:length(list)   
%            [tdislike_vs_dontcare] = kh_figure(list(k).name, filter, 'tdislike_vs_tdontcare.mat')
%             if ~isempty(tdislike_vs_dontcare)
%                tdislike_vs_dontcare_all.(list(k).name).accuracy=tdislike_vs_dontcare.tdislike_vs_tdontcare.MEG.Accuracy;
%                tdislike_vs_dontcare_all.(list(k).name).binomial=tdislike_vs_dontcare.tdislike_vs_tdontcare.MEG.Binominal;
%                tdislike_vs_dontcare_all.(list(k).name).binomial_sig = find (tdislike_vs_dontcare.tdislike_vs_tdontcare.MEG.Binominal<=0.05);
%                tdislike_vs_dontcare_all.(list(k).name).name = list(k).name;
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

end