

function adi_artifact_cleaningMEG(PathFranzi, pathInput, pathOutput, subject, filter)

%%  entfernt Kanäle, die in Datei "art_rej.mat" gespeichert sind (erstellt mittels "Auslesen_Der_artifakte.m")

FileList = dir(strcat(pathInput, subject, filesep, 'MEG_EEG_Input\noisereduced\', filter, filesep, '02_Export_Bst2FT\'));
FileList(1:2) = [];

% lade Liste mit entfernten trials/channels
load ('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\artifacts\art_rej_all_subj.mat');

for j = 1:length (FileList) % 
    load (strcat(pathInput, subject, filesep, 'MEG_EEG_Input\noisereduced\', filter, filesep, '02_Export_Bst2FT\', FileList(j).name));
    subject
    FileList(j).name 
    % figure before cleaning:
    cfgn                = [];
    cfgn.parameter      = 'trial';
    cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
    cfgn.vartrllength   = 2;
    tRetVal            = ft_timelockanalysis(cfgn, RetVal);  
    figure
    plot(tRetVal.time, tRetVal.avg(1:248,:))     
    axis tight;
    title(strcat('MEG_', FileList(j).name(1:end-4)))
    PathFigure_beforeCleaning = strcat(pathOutput, 'figures\before_cleaning\');
    
    if ~exist (PathFigure_beforeCleaning)
       mkdir(PathFigure_beforeCleaning)
    end 
    
    saveas (gcf, (strcat (PathFigure_beforeCleaning, FileList(j).name(1:end-4))), 'fig') % save MEG after cleaning figure to directory 
    close 

    
%%  entfernte Trials aus Brainstorm-Datei:
    Run = strcat('Run', FileList(j).name(end-4));
    load ( strcat(PathFranzi, subject, filesep, FileList(j).name(1:end-4), filesep, 'brainstormstudy.mat'));
   
     % new sampleInfo: 3052
    TimeSamples_new = 3052; 
    sampleinfo = zeros(length(RetVal.trial), 2);
    sampleinfo(1,:) = [1, TimeSamples_new];
    
    for k = 2:length(RetVal.trial)
        sampleinfo(k, :) = [1+TimeSamples_new, 3052*k];
        TimeSamples_new = TimeSamples_new+3052;
    end
    sampleinfo_orig = sampleinfo;
    
    if ~isempty (BadTrials) 
        trial2rej = zeros(1, length(BadTrials));
        for k =1:length(BadTrials)
            if contains(BadTrials, '__COND')
               ind_trial = cell2num(strfind(BadTrials(k), 'trial'));
               trial2rej(k) = str2double(BadTrials{k}((ind_trial+5):(ind_trial+7)));
            else
                trial2rej(k) =  str2double(BadTrials{k}((end-6):(end-4)));
            end 
        end
     
        for k = 1:length(trial2rej)
            RetVal.trial{trial2rej(k)} = [];
            RetVal.time{trial2rej(k)} = [];
            RetVal.ChannelFlag_Bst{trial2rej(k)} = [];
        end
        RetVal.trial(cellfun('isempty',RetVal.trial)) = [];
        RetVal.time(cellfun('isempty',RetVal.time)) = [];
        sampleinfo([trial2rej],:) = [];
        RetVal.ChannelFlag_Bst(cellfun('isempty',RetVal.ChannelFlag_Bst)) = [];
    end

%% entferne schlechte Kanäle (aus art_rej.mat):

    for k = 1:length(RetVal.ChannelFlag_Bst)
        BadChans = find (RetVal.ChannelFlag_Bst{1, k}(1:248) == -1);
        RetVal.trial{1,k}(BadChans, :) = NaN;
        clear BadChans;
        if  isfield (art_rej.(subject),(FileList(j).name(1:end-9)))
            if  isfield (art_rej.(subject).(FileList(j).name(1:end-9)).(Run), 'rej_channels');
                RetVal.trial{1,k}((art_rej.(subject).(FileList(j).name(1:end-9)).(Run).rej_channels),:)= NaN;
            end
        end
    end
    
     % nur MEG behalten:
     RetVal.label(272:end) = [];
     RetVal = rmfield(RetVal, 'elec');
     for k = 1:length(RetVal.trial)
         RetVal.trial{1,k}(272:end,:) = [];
     end
     
     if exist('trial2rej', 'var')
        RetVal.rejectedTrials = trial2rej;
     end
     RetVal.sampleinfo = sampleinfo ;
     RetVal.sampleinfo_orig = sampleinfo_orig;

   %% avg figure:
    
    tRetVal            = ft_timelockanalysis(cfgn, RetVal);  
    figure
    plot(tRetVal.time, tRetVal.avg(1:248,:))     
    axis tight;
    title(strcat('MEG_', FileList(j).name(1:end-4)))
    PathFigure = strcat(pathOutput, 'figures\after_cleaning\');
    if ~exist (PathFigure)
       mkdir(PathFigure)
    end 
    saveas (gcf, (strcat (PathFigure, FileList(j).name(1:end-4))), 'fig') % save MEG after cleaning figure to directory 
    close
    
   %% save clean MEG:
   cleanMEG = RetVal;
   save (strcat(pathOutput, FileList(j).name), 'cleanMEG') ;
   % clear all variables exept the initial ones:
   clearvars -except PathFranzi pathInput pathOutput subject filter FileList art_rej
   
end