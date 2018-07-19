

function adi_rejectvisual_MEG_extra
% Nachbereinigung

    subject = 'nl_adi_17'; % Probandennr. eingeben! Ansonsten einfach auf Run drücken
    filter = '1_95Hz';
    subj_mainpath       = strcat('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\', subject, '\MEG_EEG_input\noisereduced\', filter, filesep );
    path_MEG  = strcat('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\', subject, '\MEG_analysis\noisereduced\', filter, '\01_clean');
   
    if ~exist (path_MEG)
        mkdir (strcat('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\', subject, filesep, 'MEG_analysis\noisereduced', filesep, filter, '\01_clean')) ;
    end
        
    PathFigure_before_cleaning = strcat(path_MEG, filesep, 'figures\before_cleaning\');
        if ~exist (PathFigure_before_cleaning)
            mkdir(path_MEG, 'figures\before_cleaning\')
        end    
        
    PathFigure_after_cleaning = strcat(path_MEG, filesep, 'figures\after_cleaning\');
        if ~exist (PathFigure_after_cleaning)
            mkdir(path_MEG, 'figures\after_cleaning\')
        end 
        
    files   =   dir(fullfile(path_MEG, '*.mat'))
    size_files = size(files);
    cfg = [];
    
    cfg.method = 'summary';%'trial'% 'summary' %, 'trial'
%     cfg.method = 'channel';
    cfg.channel = 'MEG'; % MEG
    cfg.keepchannel = 'nan';
    cfg.latency = [-0.5 0];
    cfg.megscale = 1;
    cfg.eegscale = 5e-8;
    cfg.interactive = 'yes';
    cfg.alim     = 1e-12; 

    
    for i = 5:(size_files(1,1))
        
        load (strcat(path_MEG, filesep, files(i).name))
        
        cfgn                = [];
        cfgn.parameter      = 'trial';
        cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
    %     cfg.channel       = {'all', '-A231', '-A248', '-TRIGGER', '-RESPONSE', '-EKG', '-EOG', '-E67', '-E68'}; % kaputte bzw. störende Kanäle rauswerfen
        cfgn.vartrllength   = 2;
        
        
        %% plot and save avg before cleaning:
             
        tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
        plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG
        axis tight
        file            = files(i).name;
        
        saveas (gcf, (strcat (PathFigure_before_cleaning, file(1:end-4))), 'fig') % save EEG after cleaning figure to directory 
        
        %% rejectvisual        
        [MEG]       = ft_rejectvisual(cfg, cleanMEG); 
        [cfg_databrowser]           = ft_databrowser(cfg,MEG); 
         [MEG2]   = ft_rejectvisual(cfg, MEG); 
         
        for m = 1:length(MEG2.trial)
             MEG2.trial{1,m}(249:271,:) = cleanMEG.trial{1,m}(249:271,:);
        end
        
         %% plot and save avg after cleaning
 
        tMEG2            = ft_timelockanalysis(cfgn,MEG2);  
        figure
        plot(tMEG2.time, tMEG2.avg(1:248,:))     
        axis tight;
        title(strcat('MEG_', file(1:end-4)))
        cleanMEG = MEG2;
        saveas (gcf, (strcat (PathFigure_after_cleaning, file(1:end-4))), 'fig') % save EEG after cleaning figure to directory 
        
        save (strcat(path_MEG, filesep, (files(i).name)), 'cleanMEG'); % save mat-file after cleaning
        clear MEG2 RetVal MEG cleanMEG
               
%         %% EEG
%          cfg.channel = 274:337; % MEG
%          [EEG, trlsel,chansel] = ft_rejectvisual(cfg, RetVal); 
%         [cfg_browser] = ft_databrowser(cfg,EEG); 
%         
%         for m=1:length(MEG.trial)
%              MEG.trial{1,m}(274:337,:) = RetVal.trial{1,m}(274:337,:)
%         end
        
        
        
%         cfg.artfctdef.reject = 'complete';
%         cfg.artfctdef.feedback = 'yes';
%         MEG_clean = ft_rejectartifact(cfg, RetVal);
%         
%         badchannel = 'A231';
%         cftb=[];
%         cfgb.channel = {'all', '-A231'};
%         [data] = ft_preprocessing(cfgb, RetVal)
        
        
      
       
       
%         title( strcat (StatsFile.name, {' '}, ConfigFile.Name, {' '}, TimeWindow.TimeWindow_string) );
%         sourceplot_sign          = strcat(Path.Statistics, '\',  ConfigFile.name, '\', StatsFile.name, '_', 'sign');
%         print('-dpng', sourceplot_sign);
%          plot_ortho               = strcat( Path.SourceAnalysis, '\', ConfigFile.name, '\', 'avg_ortho', '_', ConfigFile.string );
%         saveas(gcf,plot_ortho,'fig') 
    
       
    end


end

