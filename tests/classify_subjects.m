
subj1 = 'nl_adi_06';
subj2 = 'nl_adi_08';
path_subj1 = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\' subj1 '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
path_subj2 = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\' subj2 '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];

files_subj1 = dir(fullfile(path_subj1, '*.mat'));

for p = 1:length(files_subj1)
    load([files_subj1(p).folder filesep files_subj1(p).name])
    for k=1:length(cleanMEG_interp.trial)
        cleanMEG_interp.run(k) = str2double(files_subj1(p).name(end-4));
        cleanMEG_interp.subj(k) = str2double(subj1(end-1:end));
        cleanMEG_interp.subj_label(k) = 1;
        if 1 == contains(files_subj1(p).name, 'Dislike')
            cleanMEG_interp.condition{k} = 'dislike';
            cleanMEG_interp.cond_label(k) = 2;
        elseif 1 == contains(files_subj1(p).name, 'Neu_Like')
            cleanMEG_interp.condition{k} = 'like';
            cleanMEG_interp.cond_label(k) = 1;
        end

    end
    
    if p > 1
         [~, cleanMEG_interp] = adi_check_fieldnames (var_subj1(p-1), cleanMEG_interp); 
    end

    var_subj1(p,:) = cleanMEG_interp;
    clear cleanMEG_interp

end

data_subj1 = cell2struct(cellfun(@horzcat,struct2cell(var_subj1(1)),struct2cell(var_subj1(2)), struct2cell(var_subj1(3)), struct2cell(var_subj1(4)), struct2cell(var_subj1(5)), struct2cell(var_subj1(6)), 'uni',0),fieldnames(var_subj1(1)),1);

%%
files_subj2 = dir(fullfile(path_subj2, '*.mat'));

for p = 1:length(files_subj2)
    load([files_subj2(p).folder filesep files_subj2(p).name])
    for k=1:length(cleanMEG_interp.trial)
        cleanMEG_interp.run(k) = str2double(files_subj2(p).name(end-4));
        cleanMEG_interp.subj(k) = str2double(subj2(end-1:end));
        cleanMEG_interp.subj_label(k) = 2;
        if 1 == contains(files_subj2(p).name, 'Dislike')
            cleanMEG_interp.condition{k} = 'dislike';
            cleanMEG_interp.cond_label(k) = 2;
        elseif 1 == contains(files_subj2(p).name, 'Neu_Like')
            cleanMEG_interp.condition{k} = 'like';
            cleanMEG_interp.cond_label(k) = 1;
        end

    end
    
    if p > 1
         [~, cleanMEG_interp] = adi_check_fieldnames (var_subj2(p-1), cleanMEG_interp); 
    end

    var_subj2(p,:) = cleanMEG_interp;
    clear cleanMEG_interp

end

data_subj2 = cell2struct(cellfun(@horzcat,struct2cell(var_subj2(1)),struct2cell(var_subj2(2)), struct2cell(var_subj2(3)), struct2cell(var_subj2(4)), struct2cell(var_subj2(5)), struct2cell(var_subj2(6)), 'uni',0),fieldnames(var_subj2(1)),1);

%%

data_all = cell2struct(cellfun(@horzcat,struct2cell(data_subj1),struct2cell(data_subj2), 'uni',0),fieldnames(data_subj1),1);
data_all.grad(2:end) = [];
data_all.cfg(2:end) = [];
data_all.fsample(2:end) = [];
data_all.label(:,2:end) = [];

cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'yes';
[data_all] = ft_resampledata(cfg, data_all);

cfg =[];
cfg.latency = [-0.5 1];
data_all = ft_selectdata(cfg, data_all);

for k = 1:length(data_all.trial)
    data_all.grad.label(249:end) = [];
    data_all.grad.chanori(249:end, :) = [];
    data_all.grad.chanpos(249:end, :) = [];
    data_all.grad.tra(249:end, :) = [];
    data_all.label(249:end) = [];
end

%% z-normalization:

for k = 1:length(data_all.trial) 

    for kk=1:size(data_all.trial{k},1)
        zscore.trial{k}(kk,:) = zeros(1,length(data_all.trial{k}(kk,:)));
        M = mean(data_all.trial{k}(kk,1:nearest(data_all.time{1,1},0)));
        STD = std(data_all.trial{k}(kk,1:nearest(data_all.time{1,1},0)));
        for ii=1:length(data_all.trial{k}(kk,:))
            zscore.trial{k}(kk,ii) = (data_all.trial{k}(kk,ii)- M)/STD;
        end
        clearvars M STD
    end
end

data_all.trial = zscore.trial;

%%

% fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_res_sel);
% fn_filename = fieldnames(filename);
% 
% diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);
% 
% for k=1:length(diff_fieldnames)
%     data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
% end
% fn_filename{end+1}='cfg';
% data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, fn_filename);


for k=1:length(data_all.trial)
    switch data_all.subj(k)
        case 6
            data_all.subj_label(k) = 1;
        case 8
            data_all.subj_label(k) = 2;
    end
end


for k=1:length(data_all.trial)
    switch data_all.condition{k}
        case 'dislike'
            data_all.like_dislike_label(k) = 2;
        case 'like'
            data_all.like_dislike_label(k) = 1;
    end
end


cfg_LDA = [] ;
cfg_LDA.method          = 'mvpa';
cfg_LDA.mvpa.classifier = 'lda'; %  'logreg' oder'lda' => multi-class Linear Discriminant Analysis (LDA)
cfg_LDA.mvpa.metric     = {'accuracy'; 'auc'};
% cfg_LDA.mvpa.param.lambda =  [0.0001 0.0006 0.0036 0.026 0.13 0.1 0.2 0.3 0.4 0.5 0.6 0.77 1 4.6416 27.82559 166.81 1000];
% cfg_LDA.mvpa.param.reg = 'shrink';%
cfg_LDA.mvpa.param.k = 5;
% cfg_LDA.mvpa.param.plot = 0;
% cfg_LDA.mvpa.param.repeat = 5;
% cfg_LDA.mvpa.repeat = 5;
cfg_LDA.balance = 'undersample';
cfg_LDA.cv              = 'kfold';
cfg_LDA.k               = 5;
cfg_LDA.repeat          = 5;
cfg_LDA.mvpa.normalise = 'none'; % 'none' or 'demean' or 'zscore'
cfg_LDA.latency     = [-0.5 1]; 
cfg_LDA.avgovertime     = 'no';
cfg_LDA.design          = data_all.subj_label';
tic
    stat_LDA = ft_timelockstatistics(cfg_LDA, data_all);
toc

mv_plot_result(stat_LDA.mvpa, stat_LDA.time)

savefig('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\checks\classify_subjects\accuracy.fig')

% wie gut kann Klassifikator like vs. dislike unterscheiden?
cfg_LDA.design          = data.like_dislike_label';
stat_LDA = ft_timelockstatistics(cfg_LDA, data);

mv_plot_result(stat_LDA.mvpa, stat_LDA.time)


%% helper functions

function [var, MEG_interp] = checkfields (var, MEG_interp)

%     if ~isfield(MEG_interp, 'dimord')
%         MEG_interp = setfield(MEG_interp, 'dimord', 'chan_time');
%     end

    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEG_interp);
    [field, row] = setdiff(fieldnames1, fieldnames2);
    var = rmfield(var, field);
    clear field row fieldnames1 fieldnames2
    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEG_interp);
    [field, row] = setdiff(fieldnames2, fieldnames1);
    MEG_interp = rmfield(MEG_interp, field);
    MEG_interp = orderfields(MEG_interp, var);
    
end


function [output_data_zscore] = ztrans_sensorspace(input_data)

output_data_zscore = input_data;


for p=1:length(input_data)
%     output_data_zscore(p).trial = [];
    for k=1:length(input_data(p).trial)
%         sensordata_all_subj_zscore(p).trial=cell(1, length(sensordata_all_subj(p).trial));
        for o=1:size(input_data(p).trial{k},1)
            output_data_zscore(p).trial{k}(o,:)=zeros(1,length(output_data_zscore(p).trial{k}(o,:)));
            M = mean(input_data(p).trial{k}(o,1:find(abs(input_data(p).time{k})==min(abs(input_data(p).time{k})))));
            STD = std(input_data(p).trial{k}(o,1:find(abs(input_data(p).time{k})==min(abs(input_data(p).time{k})))));
            for ii=1:length(input_data(p).trial{k}(o,:))
                output_data_zscore(p).trial{k}(o,ii) = (input_data(p).trial{k}(o,ii)- M)/STD;
            end
            clearvars M STD
        end
    end
end

clear input_data

end