function adi_select_files (inpath, subject, outpath)

    ListDir             = dir(fullfile([inpath subject filesep 'N*'])); % get a list of folders in directory
    
    for k=1:length(ListDir)
        fileList = dir(fullfile([ListDir(k).folder filesep ListDir(k).name filesep '*_average_*.mat'])); 
        if 1 == ~isempty(fileList)
            for kk = 1:length(fileList)
                delete ([fileList(kk).folder filesep fileList(kk).name])
            end
        end
    end
end



