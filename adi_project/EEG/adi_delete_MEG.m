function [] = adi_delete_MEG(inPath, outPath)

file_list = dir(fullfile(inPath, '*.mat'));
 for k = 1:length(file_list)
     if ~exist([outPath, file_list(k).name], 'file')
        load([inPath, file_list(k).name])
        for m = 1:length(RetVal.trial)
            RetVal.trial{1,m}([1:273 279 288 313 322],:) = [];
            RetVal.ChannelFlag_Bst{1,m}([1:273 279 288 313 322],:) = [];
        end
        RetVal.label ([1:273 279 288 313 322]) = [];
        RetVal = rmfield(RetVal, 'grad'); 
        EEG = RetVal;
        clear RetVal
        save ([outPath, file_list(k).name(1:end-4)], 'EEG')
        clear EEG
     end
 end
     
end