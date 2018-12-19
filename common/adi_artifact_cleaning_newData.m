

function adi_artifact_cleaningMEG(path2retval, path2cleanfile, subject, filter, trigger)

%%  entfernt Kanäle, die in Datei "art_rej.mat" gespeichert sind (erstellt mittels "Auslesen_Der_artifakte.m")

FileList = dir(fullfile([path2retval 'Neu*.mat']));

for j = 4:length (FileList) % 
    load ([FileList(j).folder filesep FileList(j).name], 'RetVal');
    
%     cfg = [];
%     cfg.trials        = 'all'; 
%     cfg.feedback = 'yes';
%     cfg.bpfilter      = 'yes'; 
%     cfg.bpfreq        = [1 45];
%     [RetVal2] = ft_preprocessing(cfg, RetVal);
%     RetVal2.ChannelFlag_Bst = RetVal.ChannelFlag_Bst;
    
    % figure before cleaning:
    cfgn                = [];
    cfgn.parameter      = 'trial';
    cfgn.baseline       = [-1 0];
    cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
    cfgn.vartrllength   = 2;
    tRetVal             = ft_timelockanalysis(cfgn, RetVal);  
    figure
    plot(tRetVal.time, tRetVal.avg(1:248,:))     
    axis tight;
    title(['MEG_', FileList(j).name(1:end-4)])
    PathFigure_beforeCleaning = [path2cleanfile 'figures\before_cleaning\'];
    
    if ~exist (PathFigure_beforeCleaning)
       mkdir(PathFigure_beforeCleaning)
    end 
    
    saveas (gcf, ([PathFigure_beforeCleaning FileList(j).name(1:end-4)]), 'fig') % save MEG after cleaning figure to directory 
    close 
%%
    cleanMEG = RetVal;
    for k = 1:length(RetVal.trial)
        ind = find(RetVal.ChannelFlag_Bst{k} == -1);
        cleanMEG.trial{k}(ind,:) = NaN; 
        cleanMEG.trial{k}(73,:) = NaN; % A231 für alle Probanden rausnehmen
        clear ind
    end

    cleanMEG = setfield(cleanMEG, 'trialinfo', []);
        
    for p = 1:length(cleanMEG.trial)
        cleanMEG.trialinfo.triggerchannel(p,:) = RetVal.trial{1,p}(272,1:end); 
        cleanMEG.trialinfo.responsechannel(p,:) = RetVal.trial{1,p}(273,1:end); 
        triggercode = find(ismember(RetVal.trial{p}(272,:),trigger.triggercodes));
        if ~isempty(triggercode)
            cleanMEG.trialinfo.triggerlabel(p) = RetVal.trial{p}(272,triggercode(end));
        elseif isempty(triggercode)
            triggercode = find(RetVal.trial{p}(272,:) >=110 & RetVal.trial{p}(272,:) ~=512);
            cleanMEG.trialinfo.triggerlabel(p) = RetVal.trial{p}(272,triggercode(end));
        end
        cleanMEG.trialinfo.response(p) = max(RetVal.trial{1,p}(273,:));
        switch cleanMEG.trialinfo.response(p)
            case 32
                cleanMEG.trialinfo.response_label{p} = 'dislike';
            case 16
                  cleanMEG.trialinfo.response_label{p} = 'like';
            case 64
                cleanMEG.trialinfo.response_label{p} = 'dontcare';
            otherwise
                cleanMEG.trialinfo.response_label{p} = FileList(j).name(1:end-9);
                switch FileList(j).name(1:end-9) 
                    case 'dislike'
                        cleanMEG.trialinfo.response(p) = 32;
                    case 'like'
                        cleanMEG.trialinfo.response(p) = 16;
                    case 'dontcare'
                        cleanMEG.trialinfo.response(p) = 64;
                end
        end
        cleanMEG.trialinfo.balldesign{p} = trigger.balldesign(find(trigger.triggercodes == cleanMEG.trialinfo.triggerlabel(p)));
        cleanMEG.trialinfo.balldesign_short{p} = trigger.balldesign_short(find(trigger.triggercodes == cleanMEG.trialinfo.triggerlabel(p)));
    end
    
    ind = find(~ismember(cleanMEG.trialinfo.triggerlabel,trigger.triggercodes));
    if ~isempty(ind)
        for p = 1:length(ind)
            Msg = ['trial nr. '  num2str(ind(p)) ' in ' FileList(j).name ' of subject ' subject ' contains triggercode no ' num2str(cleanMEG.trialinfo.triggerlabel(ind(p)))  ];
            fid = fopen(fullfile('W:\neurochirurgie\science\Kirsten\adidas\', 'wrong_triggercodes.txt'), 'a');
%             if fid == -1
%               error('Cannot open log file.');
%             end
            fprintf(fid, '%s: %s\n', datestr(now, 0), Msg);
            fclose(fid);
        end
        cleanMEG.trialinfo.triggerlabel(ind) = [];
        cleanMEG.trialinfo.balldesign(ind) = [];
        cleanMEG.trialinfo.balldesign_short(ind) = [];
        cleanMEG.trial(ind)=[];
        cleanMEG.time(ind)=[];
        cleanMEG.ChannelFlag_Bst(ind)=[];
        cleanMEG.trialinfo.triggerchannel(ind,:)=[];
        cleanMEG.trialinfo.responsechannel(ind,:)=[];
        cleanMEG.trialinfo.response(ind)=[];
        cleanMEG.trialinfo.response_label(ind)=[];
    end
    
      % nur MEG behalten:
     cleanMEG.label(249:end) = [];
     cleanMEG = rmfield(cleanMEG, 'elec');
     for k = 1:length(cleanMEG.trial)
         cleanMEG.trial{k}(249:end,:) = [];
     end

    trials = cell(1, length(cleanMEG.trial)); 
    for k=1:length(cleanMEG.trial)
        time_zero = find(cleanMEG.time{k} == 0); 
        mean_prestim = mean(cleanMEG.trial{k}(:, (time_zero/2):time_zero-5),2);
        trials{k} = cleanMEG.trial{k}-mean_prestim;
    end
    cleanMEG = rmfield(cleanMEG, 'trial');
    cleanMEG.trial = trials;

     
    %% avg figure:
    tcleanMEG = ft_timelockanalysis(cfgn, cleanMEG);  
    figure
    plot(tcleanMEG.time, tcleanMEG.avg)     
    axis tight;
    title(['MEG_', FileList(j).name(1:end-4)]);
    PathFigure = ([path2cleanfile, 'figures\after_cleaning\']);
    if ~exist (PathFigure)
       mkdir(PathFigure)
    end 
    saveas (gcf, ([PathFigure, FileList(j).name(1:end-4)]), 'fig') % save MEG after cleaning figure to directory 
    close
    
  
    save ([path2cleanfile FileList(j).name], 'cleanMEG' )
    clear cleanMEG RetVal   
end

end

   
