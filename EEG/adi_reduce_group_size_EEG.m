function adi_reduce_group_size_MEG(inPath, outPath)

%% dislike:
list_files_dislike = dir(fullfile([inPath, 'dislike*.mat']));

for i = 1:length(list_files_dislike)
        if exist([inPath, list_files_dislike(i).name])
            load([inPath, list_files_dislike(i).name])
            dislike_all_subj_red = struct('label', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', [], 'fsample', [], 'dimord',[]);

            cfg = [];
            cfg.latency = [-0.5 1];
            for k = 1:length(dislike_all_subj)
                data_out = ft_selectdata(cfg, dislike_all_subj(k))
                if ~isfield(data_out, 'sampleinfo');
                    data_out = setfield(data_out, 'sampleinfo', []);
                end
                data_out = orderfields(data_out, dislike_all_subj_red);
                dislike_all_subj_red(k) = data_out;
                clear data_out
            end
            clear dislike_all_subj
            dislike_all_subj_red(2).trial(47:66) = [];
            dislike_all_subj_red(2).time(47:66) = [];
            dislike_all_subj_red(2).sampleinfo(47:66, :) = [];
            save ([outPath, list_files_dislike(i).name(1:end-4)], 'dislike_all_subj_red', '-v7.3')
            delete([inPath, list_files_dislike(i).name])
        end
end

%% like:
list_files_like = dir(fullfile([inPath, 'like*.mat']));
like_all_subj_red = struct('label', [], 'sampleinfo', [], 'trial', [], 'time', [], 'cfg', [], 'fsample', [], 'dimord',[]);
for i = 1:length(list_files_like)
    if exist([inPath, list_files_like(i).name])
        load([inPath, list_files_like(i).name])
        % 
        cfg = [];
        cfg.latency = [-0.5 1];
        for k = 1:length(like_all_subj)
            data_out = ft_selectdata(cfg, like_all_subj(k))
            if ~isfield(data_out, 'sampleinfo');
                data_out = setfield(data_out, 'sampleinfo', []);
            end
            data_out = orderfields(data_out, like_all_subj_red);
            like_all_subj_red(k) = data_out;
            clear data_out
        end
        clear like_all_subj
        like_all_subj_red(2).trial(45:61) = [];
        like_all_subj_red(2).time(45:61) = [];
        like_all_subj_red(2).sampleinfo(45:61, :) = [];
        save ([outPath, list_files_like(i).name(1:end-4)], 'like_all_subj_red', '-v7.3')
        delete([inPath, list_files_like(i).name])
    end

end



end