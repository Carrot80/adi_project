function [data] = kh_trial2dat(trials)

if 1 == size(trials,1)

    data = zeros(length(trials), size(trials{1,1},1), size(trials{1,1},2));
    for kk = 1:length(trials)
        data(kk,:,:) = trials{kk}; 
    end
    
else
    for kk = 1:size(trials,1)
        data{kk} = squeeze(trials(kk,:,:)); 
        
    end
end

end