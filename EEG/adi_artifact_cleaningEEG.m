

function adi_artifact_cleaningEEG(PathFranzi, pathInput, pathOutput, subject)

file_list = dir(fullfile([pathInput, filesep, '*like*.mat'])); 

for i = 1:length(file_list)
    if ~exist([pathOutput, filesep, file_list(i).name])
        load([pathInput, filesep, file_list(i).name]);
        [sampleinfo] = adi_figure(EEG, pathOutput,  file_list(i).name, 'before_cleaning')

        % lade Liste mit entfernten trials/
        if contains(pathInput, 'nl_adi_11') && strncmp(file_list(i).name, 'like500_', 8)
            load ([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG_-problem_run2_and_3\brainstormstudy.mat']);
        elseif contains(pathInput, 'nl_adi_12') 
            load ([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG_problem_run123\brainstormstudy.mat']);        
        else
            load ([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\brainstormstudy.mat']);
        end
            % finde entfernte Trials
        if ~isempty (BadTrials)
           [badTrial_num, bp_num, bp_02_num] = adi_bad_trials (BadTrials, file_list(i).name(1:end-4))
            for k = 1:length(badTrial_num)
                EEG.trial{badTrial_num(k)} = [];
                EEG.time{badTrial_num(k)} = [];
                EEG.ChannelFlag_Bst{badTrial_num(k)} = []; % Channelflag enthält keine entfernten Kanäle, da ich keine Kanäle in Brainstorm-Datenbank entfernt habe; siehe Channelflag von Franzis Datei
            end

            for k = 1:length(bp_num)
                EEG.trial{bp_num(k)} = [];
                EEG.time{bp_num(k)} = [];
                EEG.ChannelFlag_Bst{bp_num(k)} = [];
            end
            for k = 1:length(bp_02_num)
                EEG.trial{bp_02_num(k)} = [];
                EEG.time{bp_02_num(k)} = [];
                EEG.ChannelFlag_Bst{bp_02_num(k)} = [];
            end

        end

        EEG.sampleinfo = sampleinfo;
        EEG.sampleinfo_orig = sampleinfo;
 
        for k = 1:length(EEG.trial)
         empty_cells(k) = isempty(EEG.trial{1, k}) 
        end

        EEG.trial(find(empty_cells == 1)) = [];
        EEG.time(empty_cells == 1) = [];
        EEG.ChannelFlag_Bst(find(empty_cells == 1)) = [];
        EEG.sampleinfo(find(empty_cells == 1), :) = [];

     

        % remove bad channels with Channelflag
        list_chans =  dir(fullfile([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\', 'data_', file_list(i).name(1:end-4), '*bandpass_02.mat']));
        if ~isempty(list_chans)
            bst_file = load([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\', list_chans(length(list_chans)).name]) 
        end
        list_chans_2 =  dir(fullfile([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\', 'data_', file_list(i).name(1:end-4), '*bandpass_bandpass.mat']));
        if ~isempty(list_chans_2) 
            bst_file = load([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\', list_chans_2(length(list_chans_2)).name]) 
        end
        list_chans_3 = dir(fullfile([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\', 'data_', file_list(i).name(1:end-4), '*_band.mat']));
        if ~isempty(list_chans_3) 
            bst_file = load([PathFranzi, subject, filesep, file_list(i).name(1:end-6), '_sensor_EEG\', list_chans_3(length(list_chans_3)).name]) 
        end
        if exist('bst_file', 'var') 
            bst_file.ChannelFlag([1:273 279 288 313 322], :) = [];

            [ind] = find(bst_file.ChannelFlag == -1);
            if ~isempty(ind)
                for k = 1:length(EEG.trial)
                    EEG.trial{1, k}(ind, :) = NaN;
                end
            end
        end
        EEG = rmfield(EEG, 'ChannelFlag_Bst');
        
        [~] = adi_figure(EEG, pathOutput, file_list(i).name, 'after_cleaning')
        
        cfg=[];
        cfg.keepchannel = 'nan';
        cfg.method = 'summary'; % oder 'trial'
        cfg.latency = [-0.5 0];
        cfg.interactive = 'yes';
        [cfg_databrowser] = ft_databrowser(cfg, EEG);
        [EEG] = ft_rejectvisual(cfg, EEG); 
        [~] = adi_figure(EEG, pathOutput, file_list(i).name, 'after_cleaning')
        [cfg_databrowser] = ft_databrowser(cfg, EEG);
        [EEG] = ft_rejectvisual(cfg, EEG); 
        
        [~] = adi_figure(EEG, pathOutput, file_list(i).name, 'after_cleaning')

       %% save clean EEG:
       cleanEEG = EEG;
       save ([pathOutput, file_list(i).name(1:end-4)], 'cleanEEG') ;
       % clear all variables exept the initial ones:
       clearvars -except PathFranzi pathInput pathOutput subject file_list i
    end
end

end

function [sampleinfo] = adi_figure(data, pathOutput, condition, cleaning_step)

cfgn                = [];
cfgn.parameter      = 'trial';
cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
cfgn.vartrllength   = 2;
tdata            = ft_timelockanalysis(cfgn, data);  
figure
plot(tdata.time, tdata.avg)     
axis tight;
title(['EEG_', condition])
PathFigure = ([pathOutput, 'figures\', cleaning_step, '\']);

if ~exist (PathFigure)
   mkdir(PathFigure)
end 

sampleinfo = tdata.sampleinfo;
saveas (gcf, ([PathFigure, condition(1:end-4)]), 'fig'); % save EEG after cleaning figure to directory 
close all

end

function [badTrial_num, bp_num, bp_02_num] = adi_bad_trials (BadTrials, condition)

 [find_BadTrials] = contains (BadTrials, condition);
 badTrials_asked_run = BadTrials(find(find_BadTrials == 1));
 % sortieren der schlechten Trials nach bandpass und ohne bandpass
 find_bp_02 = contains(badTrials_asked_run, 'bandpass_02')
 bp_02 = badTrials_asked_run(find(find_bp_02 == 1));
 badTrials_asked_run(find(find_bp_02 == 1)) = {'deleted'};
 find_bp = contains(badTrials_asked_run, 'bandpass')
 bp = badTrials_asked_run(find(find_bp == 1));
 badTrials_asked_run(find(find_bp == 1)) = {'deleted'};
 badTrials_asked_run(find(contains(badTrials_asked_run, 'deleted'))) = [];

if isempty(badTrials_asked_run)
    badTrial_num = [];
else
     for k = 1:length(badTrials_asked_run)
         if contains(badTrials_asked_run{1,k}, '__CONDdislike500')
             badTrial_num(k) = str2num(badTrials_asked_run{1,k}(end-28:end-26))
         elseif contains(badTrials_asked_run{1,k}, '__CONDlike500')
              badTrial_num(k) = str2num(badTrials_asked_run{1,k}(end-25:end-23))
         else
             badTrial_num(k) = str2num(badTrials_asked_run{1,k}(end-6:end-4))
         end
     end
end

if isempty(bp)
   bp_num = [];
else
   for k = 1:length(bp)
     if contains(badTrials_asked_run{1,k}, '__CONDdislike500')
          bp_num(k) = str2num(bp{1,k}(24:26))
     elseif contains(badTrials_asked_run{1,k}, '__CONDlike500')
          bp_num(k) = str2num(bp{1,k}(21:23))
     else
          bp_num(k) = str2num(bp{1,k}(end-15:end-13))
     end
  end
end

if isempty(bp_02)
   bp_02_num = [];
else
    for k = 1:length(bp_02)
         if contains(badTrials_asked_run{1,k}, '__CONDdislike500')
              bp_02_num(k) = str2num(bp_02{1,k}(24:26));
         elseif contains(badTrials_asked_run{1,k}, '__CONDlike500')
              bp_02_num(k) = str2num(bp_02{1,k}(21:23))
         else
              bp_02_num(k) = str2num(bp_02{1,k}(end-18:end-16))
        end
    end

end
end
