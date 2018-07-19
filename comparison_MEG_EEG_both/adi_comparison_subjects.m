function [table80] = adi_figure_comparison(table80, fieldtripPath, outpath, freq)

ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_45Hz';

if 1 == strcmp(freq, 'bp1-45Hz')
    table80.MEG.bp1_45Hz=[];
    table80.EEG.bp1_45Hz=[];
    table80.MEG_EEG.bp1_45Hz=[];
else
    table80.MEG.(freq)=[];
    table80.EEG.(freq)=[];
    table80.MEG_EEG.(freq)=[];
end


for i = 1:length(ListSubj)

    EEGpath = ([fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\04_timelockstatistics_interp\']);
    MEGpath = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\04_timelockstatistics_interp\']);
    MEG_EEG_path = ([fieldtripPath ListSubj(i).name '\combined_MEG_EEG_analysis\02c_statistics_allRuns\']);    


    EEGstats = load ([EEGpath 'like_vs_dislike_allRuns_' freq '.mat']);
    MEGstats = load ([MEGpath 'like_vs_dislike_allRuns_' freq '.mat']);
    MEG_EEG_stats = load ([MEG_EEG_path 'like_vs_dislike_allRuns_' freq '.mat']);
    
    EEGstats_subjects(i,:) =  EEGstats.Condition1vs2.Accuracy;
    MEGstats_subjects(i,:) =  MEGstats.Condition1vs2.Accuracy;
    MEG_EEG_stats_subjects(i,:) =  MEG_EEG_stats.Condition1vs2.Accuracy;
    time = MEG_EEG_stats.Condition1vs2.latency(1,:);
    clear EEGstats MEGstats MEG_EEG_stats

end

ind80 = nan(length(ListSubj), length(EEGstats_subjects));

for k=1:length(ListSubj)
   [ind] = find(MEG_EEG_stats_subjects(k,1:99) >= 0.8);
    ind80(k,ind) = MEG_EEG_stats_subjects(k,ind);
end

% figure
% plot(time(50:end), ind80(:, 50:end))
if 1==strcmp(freq, 'bp1-45Hz')
    table80.MEG_EEG.bp1_45Hz = ind80; 
else
    table80.MEG_EEG.(freq) = ind80; 
end

% length_ind = sum(isnan(table80.MEG_EEG.(freq) (:, 50:end))');
% length_intervall = length(table80.MEG_EEG.(freq))-length_ind;

% if 1==strcmp(freq, 'bp1-45Hz')
%     table80.MEG_EEG.bp1_45Hz.length_intervall = length_intervall; 
% else
% %     table80.MEG_EEG.(freq).length = length_intervall; 
% end


%%
ind_EEG_80 = nan(length(ListSubj), length(EEGstats_subjects));

for k=1:length(ListSubj)
   [ind] = find(EEGstats_subjects(k,1:99) >= 0.8);
    ind_EEG_80(k,ind) = EEGstats_subjects(k,ind);
end
figure
plot(time(50:end), ind_EEG_80(:, 50:end));

if 1==strcmp(freq, 'bp1-45Hz')
    table80.EEG.bp1_45Hz = ind_EEG_80; 
else
    table80.EEG.(freq) = ind_EEG_80; 
end

% length_ind = sum(isnan(table80.EEG.(freq) (:, 50:end))');
% length_intervall = length(table80.EEG.(freq))-length_ind;
% 
% if 1==strcmp(freq, 'bp1-45Hz')
%     table80.EEG.bp1_45Hz.length_intervall = length_intervall; 
% else
%     table80.EEG.(freq).length_intervall = length_intervall; 
% end

%% MEG
ind_MEG_80 = nan(length(ListSubj), length(MEGstats_subjects));

for k=1:length(ListSubj)
   [ind] = find(MEGstats_subjects(k,1:99) >= 0.8);
    ind_MEG_80(k,ind) = MEGstats_subjects(k,ind);
end
figure
plot(time(50:end), ind_MEG_80(:, 50:end))

if 1==strcmp(freq, 'bp1-45Hz')
    table80.MEG.bp1_45Hz = ind_MEG_80; 
else
    table80.MEG.(freq) = ind_MEG_80; 
end
% 
% length_ind = sum(isnan(table80.MEG.(freq) (:, 50:end))');
% length_intervall = length(table80.MEG.(freq))-length_ind;

% if 1==strcmp(freq, 'bp1-45Hz')
%     table80.MEG.bp1_45Hz.length_intervall = length_intervall; 
% else
%     table80.MEG.(freq).length_intervall = length_intervall; 
% end

end