
function  adi_appenddata_freq(inPath, outPath, freqname)

% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
         
% dislike
if ~exist(strcat(outPath, 'dislike_allRuns_', freqname, '.mat'), 'file');
    fileName = strcat('*dislike*', freqname, '.mat');
    files_dislike = dir(fullfile(inPath, fileName));
    size_files_dislike = size(files_dislike);

    for i = 1:(size_files_dislike(1,1))
        load (strcat(inPath, files_dislike(i).name));
        files2append(:,i) = data_bpfreq;
        clear data_bpfreq
    end

    cfg = [];
    switch size_files_dislike(1)
        case 3
            dislike_allRuns = ft_appenddata(cfg, files2append(1), files2append(2), files2append(3))
            % dislike_allRuns=setfield(dislike_allRuns, 'grad', files2append(1).grad)
        case 2
             dislike_allRuns = ft_appenddata(cfg, files2append(1), files2append(2))
        %      dislike_allRuns=setfield(dislike_allRuns, 'grad', files2append(1).grad)
    end
     save (strcat(outPath, 'dislike_allRuns_', freqname, '.mat'), 'dislike_allRuns');
     clear files2append
end

%% like 
if ~exist(strcat(outPath, filesep, 'like_allRuns_', freqname, '.mat'));
    fileName = strcat('like*', freqname, '.mat');
    files_like = dir(fullfile(inPath, fileName));
    size_files_like = size(files_like);

    for i = 1:(size_files_like(1,1))
        load (strcat(inPath, files_like(i).name));
        files2append(:,i)= data_bpfreq;
        clear data_bpfreq;    
    end

    cfg=[];
    switch size_files_like(1) 
        case 3
        like_allRuns = ft_appenddata(cfg, files2append(1), files2append(2), files2append(3))
    %     like_allRuns=setfield(like_allRuns, 'grad', files2append(1).grad)

        case 2
        like_allRuns = ft_appenddata(cfg, files2append(1), files2append(2))
    %     like_allRuns=setfield(like_allRuns, 'grad', files2append(1).grad)
    end

    save (strcat(outPath, filesep, 'like_allRuns_', freqname, '.mat'), 'like_allRuns');
    clear files2append 
end

  %% dontcare
 if ~exist(strcat(outPath, filesep, 'dontcare_allRuns_', freqname, '.mat'));
    fileName = strcat('*dontcare*', freqname, '.mat');
    files_dontcare = dir(fullfile(inPath, fileName)); 
    size_files_dontcare = size(files_dontcare);

    if size_files_dontcare(1,1) == 0
        return
    else
        for i=1:size_files_dontcare(1,1)
            load(strcat(inPath, files_dontcare(i).name))
%              if ~isfield(data_bpfreq, 'dimord')
%                     data_bpfreq = setfield(data_bpfreq, 'dimord', 'chan_time');
%                     data_bpfreq = orderfields(data_bpfreq, files2append(1));
%              end
             files2append(:,i)= data_bpfreq;
             clear data_bpfreq;      
        end

        cfg=[];
        switch size_files_dontcare(1) 
            case 3
                dontcare_allRuns = ft_appenddata(cfg, files2append(1), files2append(2), files2append(3))
        %       dontcare_allRuns=setfield(dontcare_allRuns, 'grad', files2append(1).grad)
                save (strcat(outPath, filesep, 'dontcare_allRuns_', freqname, '.mat'), 'dontcare_allRuns');

            case 2
                dontcare_allRuns = ft_appenddata(cfg, files2append(1), files2append(2))
        %       dontcare_allRuns=setfield(dontcare_allRuns, 'grad', files2append(1).grad)
                save (strcat(outPath, filesep, 'dontcare_allRuns_', freqname, '.mat'), 'dontcare_allRuns');
            case 1
                dontcare_allRuns = files2append; 
                save (strcat(outPath, filesep, 'dontcare_allRuns_', freqname, '.mat'), 'dontcare_allRuns');
        end
    end
  end
end