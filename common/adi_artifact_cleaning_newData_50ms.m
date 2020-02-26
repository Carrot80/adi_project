

function adi_artifact_cleaningMEG(path2ft_data, bst_path, path2cleanfile, subject, filter, trigger)

%% Achtung: benutzt nicht infos, die in Datei "art_rej.mat" gespeichert sind (erstellt mittels "Auslesen_Der_artifakte.m")

FileList = dir(fullfile([path2ft_data 'Neu*ike50_*.mat']));

for jj = 1:length (FileList) % 
    if ~exist ([path2cleanfile FileList(jj).name], 'file')
        load ([FileList(jj).folder filesep FileList(jj).name], 'ft_data');
        load ([bst_path FileList(jj).name(1:end-4) filesep 'brainstormstudy.mat'])
               
        if ~isempty(BadTrials)
            for ii = 1:length(BadTrials)
                bad(ii) = str2num(BadTrials{ii,1}(end-6:end-4));
            end
            ft_data.trial(bad) = [];
            ft_data.time(bad) = [];
            ft_data.ChannelFlag_Bst(bad) = [];
            clear bad
        end
        
        cleanMEG = ft_data;
        for kk = 1:length(ft_data.trial)
            ind = find(ft_data.ChannelFlag_Bst{kk} == -1);
            cleanMEG.trial{kk}(ind,:) = NaN; 
            clear ind
        end
        label_A231 = find(strcmp(cleanMEG.label, 'A231')); % A231 für alle Probanden rausnehmen
        
        for kk = 1:length(cleanMEG.trial)
            cleanMEG.trial{kk}(label_A231,:) = NaN; 
        end
        
        % figure before cleaning:
%         cfgn                = [];
%         cfgn.parameter      = 'trial';
%         cfgn.baseline       = [-0.5 0];
%         cfgn.vartrllength   = 2;
%         t_cleanMEG             = ft_timelockanalysis(cfgn, cleanMEG);  
%         figure
%         plot(cleanMEG.time{1,1}, t_cleanMEG.avg(1:248,:))     
%         axis tight;
%         title(['MEG_', FileList(jj).name(1:end-4)])
%         PathFigure_beforeCleaning = [path2cleanfile 'figures\before_cleaning\'];
% 
%         if ~exist (PathFigure_beforeCleaning)
%            mkdir(PathFigure_beforeCleaning)
%         end 
% 
%         saveas (gcf, ([PathFigure_beforeCleaning FileList(jj).name(1:end-4)]), 'fig') % save MEG after cleaning figure to directory 
%         close 
    
   
        
        %%
        cleanMEG = setfield(cleanMEG, 'trialinfo', []);

        for p = 1:length(cleanMEG.trial)
            cleanMEG.trialinfo.triggerchannel(p,:) = ft_data.trial{1,p}(272,1:end); 
            cleanMEG.trialinfo.responsechannel(p,:) = ft_data.trial{1,p}(273,1:end); 
            triggercode = find(ismember(ft_data.trial{p}(272,:),trigger.triggercodes));
            if ~isempty(triggercode)
                cleanMEG.trialinfo.triggerlabel(p) = ft_data.trial{p}(272,triggercode(end));
            elseif isempty(triggercode)
                triggercode = find(ft_data.trial{p}(272,:) >=110 & ft_data.trial{p}(272,:) ~=512);
                cleanMEG.trialinfo.triggerlabel(p) = ft_data.trial{p}(272,triggercode(end));
            end
            cleanMEG.trialinfo.response(p) = max(ft_data.trial{1,p}(273,:));
            switch cleanMEG.trialinfo.response(p)
                case 32
                    cleanMEG.trialinfo.response_label{p} = 'dislike';
                case 16
                      cleanMEG.trialinfo.response_label{p} = 'like';
                case 64
                    cleanMEG.trialinfo.response_label{p} = 'dontcare';
                otherwise
                    cleanMEG.trialinfo.response_label{p} = FileList(jj).name(1:end-9);
                    switch FileList(jj).name(1:end-9) 
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
                Msg = ['trial nr. '  num2str(ind(p)) ' in ' FileList(jj).name ' of subject ' subject ' contains triggercode no ' num2str(cleanMEG.trialinfo.triggerlabel(ind(p)))  ];
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
        

        %% avg figure:
        
        cfgn                = [];
        cfgn.parameter      = 'trial';
        cfgn.baseline       = [-0.5 0];
        cfgn.vartrllength   = 2;
        tcleanMEG = ft_timelockanalysis(cfgn, cleanMEG);  
        figure
        plot(tcleanMEG.time, tcleanMEG.avg)     
        axis tight;
        title(['MEG_', FileList(jj).name(1:end-4)]);
        PathFigure = ([path2cleanfile, 'figures\after_cleaning\']);
        if ~exist (PathFigure)
           mkdir(PathFigure)
        end 
        saveas (gcf, ([PathFigure, FileList(jj).name(1:end-4)]), 'fig') % save MEG after cleaning figure to directory 
        close


        save ([path2cleanfile FileList(jj).name], 'cleanMEG' )
        clear cleanMEG ft_data   
     end
end
end

   
