

function adi_rejectvisual_MEG_extra(path2cleanfile, subject)

% Nachbereinigung
    PathFigure_after_cleaning = [path2cleanfile, 'figures\after_cleaning\'];   

    files   =   dir(fullfile(path2cleanfile, '*.mat'));
    size_files = size(files);
    
    for ii = 1:(size_files(1,1))
        
        load ([files(ii).folder, filesep, files(ii).name], 'cleanMEG')
        
        if ~isfield(cleanMEG, 'additional_cleaning')
            cfgn                = [];
            cfgn.parameter      = 'trial';
            cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
            cfgn.vartrllength   = 2;


            %% plot and save avg before cleaning:

            tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
            figure
            plot(tcleanMEG.time, squeeze(mean(tcleanMEG.trial(:,1:248,:),1))) % MEG
            axis tight    

            %% rejectvisual       

%             cfg = [];
%             out_databrowser = ft_databrowser(cfg, cleanMEG);

             if 1 == strcmp(subject, 'nl_adi_19') && ~isfield(cleanMEG, 'additional_cleaning') 
                bad_chan = {'A248'; 'A246'};
                bad_chan_ind = [];
                for kk = 1:length(bad_chan)
                     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
                end

                for kk = 1:length(cleanMEG.trial)
                    cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
                end
                
             
            end
             if 1 == strcmp(subject, 'nl_adi_20') && ~isfield(cleanMEG, 'additional_cleaning') 
                 
                bad_chan = {'A248'; 'A246'; 'A247'; 'A228'};
                bad_chan_ind = [];
                for kk = 1:length(bad_chan)
                     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
                end

                for kk = 1:length(cleanMEG.trial)
                    cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
                end
               
                
             end

            if 1 == strcmp(subject, 'nl_adi_21') && ~isfield(cleanMEG, 'additional_cleaning') % entferne A22, A248 und A247
                bad_chan = {'A22'; 'A247'; 'A248'; 'A231'; 'A228'};
                bad_chan_ind = [];
                for kk = 1:length(bad_chan)
                     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
                end

                for kk = 1:length(cleanMEG.trial)
                    cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
                end
                if 1 == strcmp(files(ii).name, 'Neu_Dislike50_1.mat')
                    cleanMEG.trial(3) = [];
                    cleanMEG.time(3) = [];
                    cleanMEG.trialinfo.triggerchannel(3,:) = [];
                    cleanMEG.trialinfo.responsechannel(3,:) = [];
                    cleanMEG.trialinfo.triggerlabel(3) = [];
                    cleanMEG.trialinfo.response(3) = [];
                    cleanMEG.trialinfo.response_label(3) = [];
                    cleanMEG.trialinfo.balldesign(3) = [];
                    cleanMEG.trialinfo.balldesign_short(3) = [];
                elseif 1 == strcmp(files(ii).name, 'Neu_Dislike50_2.mat') || 1 == strcmp(files(ii).name, 'Neu_Dislike50_3.mat')
                    cleanMEG.trial(11) = [];
                    cleanMEG.time(11) = [];
                    cleanMEG.trialinfo.triggerchannel(11,:) = [];
                    cleanMEG.trialinfo.responsechannel(11,:) = [];
                    cleanMEG.trialinfo.triggerlabel(11) = [];
                    cleanMEG.trialinfo.response(11) = [];
                    cleanMEG.trialinfo.response_label(11) = [];
                    cleanMEG.trialinfo.balldesign(11) = [];
                    cleanMEG.trialinfo.balldesign_short(11) = [];
                elseif 1 == strcmp(files(ii).name, 'Neu_Like50_1.mat')
                    cleanMEG.trial([1, 2]) = [];
                    cleanMEG.time([1, 2]) = [];
                    cleanMEG.trialinfo.triggerchannel([1, 2],:) = [];
                    cleanMEG.trialinfo.responsechannel([1, 2],:) = [];
                    cleanMEG.trialinfo.triggerlabel([1, 2]) = [];
                    cleanMEG.trialinfo.response([1, 2]) = [] ;
                    cleanMEG.trialinfo.response_label([1,2]) = [];
                    cleanMEG.trialinfo.balldesign([1, 2]) = [];
                    cleanMEG.trialinfo.balldesign_short([1, 2]) = [];
                 elseif 1 == strcmp(files(ii).name, 'Neu_Like50_2.mat')
                    cleanMEG.trial(33) = [];
                    cleanMEG.time(33) = [];
                    cleanMEG.trialinfo.triggerchannel(33,:) = [];
                    cleanMEG.trialinfo.responsechannel(33,:) = [];
                    cleanMEG.trialinfo.triggerlabel(33) = [];
                    cleanMEG.trialinfo.response(33) = [];
                    cleanMEG.trialinfo.response_label(33) = [];
                    cleanMEG.trialinfo.balldesign(33) = [];
                    cleanMEG.trialinfo.balldesign_short(33) = [];
                    
                elseif 1 == strcmp(files(ii).name, 'Neu_Like50_2.mat')
                    cleanMEG.trial([3, 13]) = [];
                    cleanMEG.time([3, 13]) = [];
                    cleanMEG.trialinfo.triggerchannel([3, 13],:) = [];
                    cleanMEG.trialinfo.responsechannel([3, 13],:) = [];
                    cleanMEG.trialinfo.triggerlabel([3, 13]) = [];
                    cleanMEG.trialinfo.response([3, 13]) = [] ;
                    cleanMEG.trialinfo.response_label([3, 13]) = [];
                    cleanMEG.trialinfo.balldesign([3, 13]) = [];
                    cleanMEG.trialinfo.balldesign_short([3, 13]) = [];
                    
                end
                
                cleanMEG.additional_cleaning = 'yes';
            end



             %% plot and save avg after cleaning

            tcleanMEG            = ft_timelockanalysis(cfgn,cleanMEG);  
            figure
            plot(tcleanMEG.time, squeeze(mean(tcleanMEG.trial(:,1:248,:),1))) % MEG
            axis tight;
            title(files(ii).name(1:end-4))

            saveas (gcf, [PathFigure_after_cleaning files(ii).name(1:end-4) '_extracleaning'], 'fig') % save EEG after cleaning figure to directory     
            save ([path2cleanfile files(ii).name], 'cleanMEG'); % save mat-file after cleaning
        end
        
        clearvars cleanMEG

    end


end

