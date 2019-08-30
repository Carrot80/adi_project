function [] = adi_ica(path2data, path2save)


%% run 1:
adi_ica(path2data, path2save, num2str(1))    


%% run2:
if ~exist([outPath_extdisc 'virtsens\run2\virtsens_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, condition, outPath_extdisc)
end

%% run3:
if ~exist([outPath_extdisc 'virtsens\run3\virtsens_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, condition, outPath_extdisc)
end

end


function [] = adi_ica(path2data, path2save, run)   

list = dir(fullfile([path2data, '*' run '.mat'])); 

for kk = 1:length(list)
    load([list(kk).folder filesep list(kk).name ]);
    if exist('RetVal', 'var')
        ft_data = RetVal;
        clear RetVal
    end
    ft_data.filename = list(kk).name;
    ft_data.cfg = [];
    data(kk) = ft_data;
    clear ft_data
end

%% vor concatenate triallänge angleichen (noch machen)
cfg = [];
cfg.latency = [-1.5 2];
for kk = 1:length(data)
    if length(data(kk).time{1}) == 4071
        data(kk) = ft_selectdata(cfg, data(kk));
    end
    if data(kk).fsample == 1017.5
        data(kk).fsample = 1017.3;
    end
end

data = setfield(data, {1}, 'condition', []);

for kk = 1:length(data)
    for pp = 1:length(data(kk).trial)
        data(kk).condition{pp} = data(kk).filename;
    end
end

data_all = data(1);
for kk = 2:length(data)
    data_all.trial = cat(2, data_all.trial, data(kk).trial);
    data_all.time = cat(2, data_all.time, data(kk).time);   
    data_all.condition = cat(2, data_all.condition, data(kk).condition); 
    data_all.label = cat(2, data_all.label, data(kk).label); 
    data_all.fsample = cat(2, data_all.fsample, data(kk).fsample); 
    data_all.ChannelFlag_Bst = cat(2, data_all.ChannelFlag_Bst, data(kk).ChannelFlag_Bst); 
end

%% sanity checks:




%%


Z = cell2struct(cellfun(@horzcat,struct2cell(data(1)),struct2cell(data(2)), struct2cell(data(3)), struct2cell(data(4)),'uni',0),fieldnames(data(1)),1);
Z = cell2struct(cellfun(@horzcat,struct2cell(data(1)),struct2cell(data(2)),'uni',0),fieldnames(data(1)),1);


%%





%% A71 und A231 interpolieren: 

% [neighbours] = MEG_neighbours (cleanMEG); 
% 
% for kk = 1:length(data)
%     ind = find(strcmp(data(kk).label, 'A71'));
%     for pp = 1:length(data(kk).trial)
%         
%         
%         
%     end
% 
% end

%%

for kk = 1:length(data)
    cfg = [];
    cfg.resamplefs = 300;
    cfg.detrend    = 'no';
    data_res(kk) = ft_resampledata(cfg, data(kk));

end

Z = cell2struct(cellfun(@vertcat,struct2cell(data_res(1)),struct2cell(data_res(2)),'uni',0),fieldnames(data_res(1)),1);
comp = ft_componentanalysis(cfg, data_res);






end




