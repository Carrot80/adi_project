function [output_data_zscore] = adi_ztrans_sensorspace(input_data)

  
output_data_zscore = input_data;

for pp=1:length(input_data)
    output_data_zscore(pp).trial = [];
    output_data_zscore(pp).trial = cell(1, length(input_data(pp).trial));
    for kk = 1:length(input_data(pp).trial)
        output_data_zscore(pp).trial{kk} = zeros(size(input_data(pp).trial{kk},1),size(input_data(pp).trial{kk},2));
        for oo = 1:size(input_data(pp).trial{kk},1)   
            M = mean(input_data(pp).trial{kk}(oo,1:find(abs(input_data(pp).time{kk}) == min(abs(input_data(pp).time{kk})))));
            STD = std(input_data(pp).trial{kk}(oo,1:find(abs(input_data(pp).time{kk}) == min(abs(input_data(pp).time{kk})))));
            output_data_zscore(pp).trial{kk}(oo,:) = (input_data(pp).trial{kk}(oo,:)- M)/STD;
            clearvars M STD
        end
    end
end

clear input_data




end







% ind_like=find(session.labels==1);
% like=session.data(ind_like,:,:);
% ind_dislike=find(session.labels==2);
% dislike=session.data(ind_dislike,:,:);
% 
% ind_zero=find(abs(session.time{1,1})==min(abs(session.time{1})));
% 
% mean_like = mean(like(:,:,1:ind_zero),3);
% mean_like_trls = mean(mean_like);
% std_like = std(like);
% std_like = std(like(1,:,1:ind_zero),3);
% mean_like_trls = mean(mean_like);
% 
% 
% for k=1:length(like)
%      M = mean(like(o,1:find(abs(input_data(p).time{k})==min(abs(input_data(p).time{k})))));
%     
% end


% for p=1:length(input_data)
% %     output_data_zscore(p).trial = [];
%     for k=1:length(input_data(p).trial)
% %         sensordata_all_subj_zscore(p).trial=cell(1, length(sensordata_all_subj(p).trial));
%         for o=1:size(input_data(p).trial{k},1)
%             output_data_zscore(p).trial{k}(o,:)=zeros(1,length(output_data_zscore(p).trial{k}(o,:)));
%             M = mean(input_data(p).trial{k}(o,1:find(abs(input_data(p).time{k})==min(abs(input_data(p).time{k})))));
%             STD = std(input_data(p).trial{k}(o,1:find(abs(input_data(p).time{k})==min(abs(input_data(p).time{k})))));
%             for ii=1:length(input_data(p).trial{k}(o,:))
%                 output_data_zscore(p).trial{k}(o,ii) = (input_data(p).trial{k}(o,ii)- M)/STD;
%             end
%             clearvars M STD
%         end
%     end
% end
% 
% clear input_data
