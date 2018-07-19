function adi_select_files_EEG (inpath, outpath, run )

    ListDir             = dir(inpath); % get a list of folders in directory
    ListDir(1:3)        = []; % remove the first three ./..  and "brainstorm.mat"  
    isTimelock          = 1; % isTimelock --> RawData (0) or TimeLock Data (1)
    cd (inpath)  
    
    ListDislike         = dir('data_dislike500*.mat'); %  nur Dislike-MatFiles ausgeben
    ChannelFile     = load ('channel_4d_acc1.mat'); %  Configuration of MEG Channels => von Brainstorm vorher exportieren
    [RetVal]        = adi_mk_ft_structure(ChannelFile, ListDislike, isTimelock) ;  % Convertiert Brainstorm-files in Fieldtrip-Files
    save ([outpath, 'dislike500_', run, '.mat'], 'RetVal');
    clear RetVal;

    ListLike         = dir('data_like500*.mat'); %  nur like-MatFiles ausgeben
    [RetVal]        = adi_mk_ft_structure(ChannelFile, ListLike, isTimelock) ;  % Convertiert Brainstorm-files in Fieldtrip-Files
    save (strcat(outpath, filesep, 'like500_', run, '.mat'), 'RetVal');
    clear RetVal;

    ListDontcare         = dir('data_dontcare500*.mat'); %  nur dontcare-MatFiles ausgeben
    if ~isempty(ListDontcare) 
        [RetVal]        = adi_mk_ft_structure(ChannelFile, ListDontcare, isTimelock) ;
        save (strcat(outpath, filesep, 'dontcare500_', run, '.mat'), 'RetVal');
        clear RetVal
    end
        
        
end