function [output_data_zscore] = adi_ztrans_sensorspace(input_data)

  
output_data_zscore = input_data;

for pp=1:length(input_data)
    output_data_zscore(pp).trial = [];
    output_data_zscore(pp).trial = cell(1, length(input_data(pp).trial));
    for kk = 1:length(input_data(pp).trial)
        output_data_zscore(pp).trial{kk} = zeros(size(input_data(pp).trial{kk},1),size(input_data(pp).trial{kk},2));
        for oo = 1:size(input_data(pp).trial{kk},1)   
            M = mean(input_data(pp).trial{kk}(oo,nearest(input_data(pp).time{kk}, -0.5):nearest(input_data(pp).time{kk}, -0.03)));
            STD = std(input_data(pp).trial{kk}(oo,nearest(input_data(pp).time{kk}, -0.5):nearest(input_data(pp).time{kk}, -0.03)));
            output_data_zscore(pp).trial{kk}(oo,:) = (input_data(pp).trial{kk}(oo,:)- M)./STD;
            clearvars M STD
        end
    end
end

clear input_data




end



% 
% 
% output_data_zscore = input_data;
% 
% for pp=1:length(input_data)
%     output_data_zscore(pp).trial = [];
%     output_data_zscore(pp).trial = cell(1, length(input_data(pp).trial));
%     for kk = 1:length(input_data(pp).trial)
%         output_data_zscore(pp).trial{kk} = zeros(size(input_data(pp).trial{kk},1),size(input_data(pp).trial{kk},2));
%         for oo = 1:size(input_data(pp).trial{kk},1)   
%             M = mean(input_data(pp).trial{kk}(oo,1:find(abs(input_data(pp).time{kk}) == min(abs(input_data(pp).time{kk})))));
%             STD = std(input_data(pp).trial{kk}(oo,1:find(abs(input_data(pp).time{kk}) == min(abs(input_data(pp).time{kk})))));
%             output_data_zscore(pp).trial{kk}(oo,:) = (input_data(pp).trial{kk}(oo,:)- M)/STD;
%             clearvars M STD
%         end
%     end
% end
% 
% clear input_data
% 
