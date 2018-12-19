

function adi_rejectvisual_MEG_extra(path2cleanfile, path_interpolated, subject)
% Nachbereinigung
    PathFigure_before_cleaning = [path2cleanfile, 'figures\before_cleaning\'];   
    PathFigure_after_cleaning = [path2cleanfile, 'figures\after_cleaning\'];   

    files   =   dir(fullfile(path2cleanfile, '*.mat'));
    size_files = size(files);
    cfg = [];
    cfg.method = 'trial';%'trial'% 'summary' %, 'trial'
%     cfg.method = 'channel';
    cfg.channel = 'MEG'; % MEG
    cfg.keepchannel = 'nan';
    cfg.latency = [-0.5 0];
    cfg.megscale = 1;
    cfg.eegscale = 5e-8;
    cfg.interactive = 'yes';
    cfg.alim     = 1e-12; 
    cfg.keeptrial = 'yes';

    
    for i = 1:(size_files(1,1))
        
        load ([files(i).folder, filesep, files(i).name], 'cleanMEG')
        load([path_interpolated files(i).name], 'RetVal')
        cfgn                = [];
        cfgn.parameter      = 'trial';
        cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
        cfgn.vartrllength   = 2;
        
        
        %% plot and save avg before cleaning:
             
        tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
        figure
        plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG
        axis tight
        file            = files(i).name;
        
%         saveas (gcf, [PathFigure_before_cleaning, file(1:end-4)], 'fig') % save EEG after cleaning figure to directory 
        
        %% rejectvisual        

        trialinfo = cleanMEG.trialinfo;
        
        if ~isequal(length(cleanMEG.trial), 1)
            [cleanMEG]       = ft_rejectvisual(cfg, cleanMEG); 
            cleanMEG.trialinfo = trialinfo;
            for p = 1:length(cleanMEG.trial)
                ind = find(isnan(cleanMEG.trial{p}(:,1)));
                cleanMEG.ChannelFlag_Bst{p}(ind) =-1;
                RetVal.ChannelFlag_Bst{p}(ind) =-1;
            end
            if ~isempty(cleanMEG.cfg.artfctdef.trial.artifact)
                for p=1:length(cleanMEG.cfg.artfctdef.trial.artifact)
                   ind(p) = find(cleanMEG.cfg.artfctdef.trial.artifact(p,1) == cleanMEG.sampleinfo(:,1));
                end
                cleanMEG.trial(find(ind))=[];
                cleanMEG.time(find(ind))=[];
                cleanMEG.sampleinfo(find(ind),:)=[];
                cleanMEG.ChannelFlag_Bst(find(ind),:)=[];
                cleanMEG.triggerchannel
                cleanMEG.
                cleanMEG
                cleanMEG
                cleanMEGcleanMEG
            end
        elseif isequal(length(cleanMEG.trial), 1)

            z_cleanMEG = zscore(cleanMEG.trial{1});
            ind=zeros(1, size(z_cleanMEG,1));
            figure
            plot(tcleanMEG.time, z_cleanMEG) % MEG
            for p=1:size(z_cleanMEG,1)
               [temp] =  find(abs(z_cleanMEG(p,1000:2000)) >=5 );
               if ~isempty(temp)
                   ind(p)=sum(temp);
               end
            end
            badchan = find(ind);
            cleanMEG.trial{1}(badchan,:)=NaN;
            cleanMEG.trialinfo = trialinfo
        end
        
        
         %% plot and save avg after cleaning
 
        tcleanMEG            = ft_timelockanalysis(cfgn,cleanMEG);  
        figure
        plot(tcleanMEG.time, tcleanMEG.avg(1:248,:))     
        axis tight;
        title(strcat('MEG_', file(1:end-4)))
        
        saveas (gcf, [PathFigure_after_cleaning file(1:end-4)], 'fig') % save EEG after cleaning figure to directory 
        
%         load ([path_interpolated, files(i).name])
        
        
        save ([path2cleanfile (files(i).name)], 'cleanMEG'); % save mat-file after cleaning
        save ([path_interpolated (files(i).name)], 'RetVal'); % save mat-file after cleaning
           clearvars cleanMEG    
  
        
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

