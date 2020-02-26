

function adi_rejectvisual_MEG_extra(path2cleanfile, subject)

% Nachbereinigung
    PathFigure_after_cleaning = [path2cleanfile, 'figures\after_cleaning\'];   

    files   =   dir(fullfile(path2cleanfile, '*.mat'));
    size_files = size(files);
    
    for ii = 8:(size_files(1,1))
        
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
         
            
             if 1 == strcmp(subject, 'nl_adi_20') && ~isfield(cleanMEG, 'additional_cleaning') 
                 
                bad_chan = {'A248'; 'A246'; 'A247'; 'A228'};
                bad_chan_ind = [];
                for kk = 1:length(bad_chan)
                     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
                end

                for kk = 1:length(cleanMEG.trial)
                    cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
                end
                cfg = [];
                out_databrowser = ft_databrowser(cfg, cleanMEG); 
                    
             if 1 == strcmp(files(ii).name, 'Neu_Dislike500_2.mat')
                    cleanMEG.trial([15, 20, 34, 35]) = [];
                    cleanMEG.time([15, 20, 34, 35]) = [];
                    cleanMEG.trialinfo.triggerchannel([15, 20, 34, 35],:) = [];
                    cleanMEG.trialinfo.responsechannel([15, 20, 34, 35],:) = [];
                    cleanMEG.trialinfo.triggerlabel([15, 20, 34, 35]) = [];
                    cleanMEG.trialinfo.response([15, 20, 34, 35]) = [];
                    cleanMEG.trialinfo.response_label([15, 20, 34, 35]) = [];
                    cleanMEG.trialinfo.balldesign([15, 20, 34, 35]) = [];
                    cleanMEG.trialinfo.balldesign_short([15, 20, 34, 35]) = [];
             elseif  1 == strcmp(files(ii).name, 'Neu_Dislike500_4.mat')
                    cleanMEG.trial([30]) = [];
                    cleanMEG.time([30]) = [];
                    cleanMEG.trialinfo.triggerchannel([30],:) = [];
                    cleanMEG.trialinfo.responsechannel([30],:) = [];
                    cleanMEG.trialinfo.triggerlabel([30]) = [];
                    cleanMEG.trialinfo.response([30]) = [];
                    cleanMEG.trialinfo.response_label([30]) = [];
                    cleanMEG.trialinfo.balldesign([30]) = [];
                    cleanMEG.trialinfo.balldesign_short([30]) = [];
              elseif  1 == strcmp(files(ii).name, 'Neu_Like500_2.mat')
                    cleanMEG.trial([11, 15, 17]) = [];
                    cleanMEG.time([11, 15, 17]) = [];
                    cleanMEG.trialinfo.triggerchannel([11, 15, 17],:) = [];
                    cleanMEG.trialinfo.responsechannel([11, 15, 17],:) = [];
                    cleanMEG.trialinfo.triggerlabel([11, 15, 17]) = [];
                    cleanMEG.trialinfo.response([11, 15, 17]) = [];
                    cleanMEG.trialinfo.response_label([11, 15, 17]) = [];
                    cleanMEG.trialinfo.balldesign([11, 15, 17]) = [];
                    cleanMEG.trialinfo.balldesign_short([11, 15, 17]) = [];
                 
             end
                
             end

            if 1 == strcmp(subject, 'nl_adi_19') && ~isfield(cleanMEG, 'additional_cleaning') 
                bad_chan = {'A248'};
                bad_chan_ind = [];
                for kk = 1:length(bad_chan)
                     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
                end

                for kk = 1:length(cleanMEG.trial)
                    cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
                end
                
             if 1 == strcmp(files(ii).name, 'Neu_Dontcare500_1.mat')
                    cleanMEG.trial(2) = [];
                    cleanMEG.time(2) = [];
                    cleanMEG.trialinfo.triggerchannel(2,:) = [];
                    cleanMEG.trialinfo.responsechannel(2,:) = [];
                    cleanMEG.trialinfo.triggerlabel(2) = [];
                    cleanMEG.trialinfo.response(2) = [];
                    cleanMEG.trialinfo.response_label(2) = [];
                    cleanMEG.trialinfo.balldesign(2) = [];
                    cleanMEG.trialinfo.balldesign_short(2) = [];
             end
            end
            
            
            if 1 == strcmp(subject, 'nl_adi_21') && ~isfield(cleanMEG, 'additional_cleaning') % entferne A22, A248 und A247
                bad_chan = {'A22'; 'A212'; 'A246'; 'A247'; 'A248'; 'A231'; 'A228'};
                bad_chan_ind = [];
                for kk = 1:length(bad_chan)
                     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
                end

                for kk = 1:length(cleanMEG.trial)
                    cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
                end
                
                cfg = [];
                out_databrowser = ft_databrowser(cfg, cleanMEG);

                
                
                if 1 == strcmp(files(ii).name, 'Neu_Dislike500_1.mat')
                    cleanMEG.trial(5) = [];
                    cleanMEG.time(5) = [];
                    cleanMEG.trialinfo.triggerchannel(5,:) = [];
                    cleanMEG.trialinfo.responsechannel(5,:) = [];
                    cleanMEG.trialinfo.triggerlabel(5) = [];
                    cleanMEG.trialinfo.response(5) = [];
                    cleanMEG.trialinfo.response_label(5) = [];
                    cleanMEG.trialinfo.balldesign(5) = [];
                    cleanMEG.trialinfo.balldesign_short(5) = [];
                elseif 1 == strcmp(files(ii).name, 'Neu_Dislike500_2.mat')
                    cleanMEG.trial([3,10,11,12,13,20,22,25,26]) = [];
                    cleanMEG.time([3,10,11,12,13,20,22,25,26]) = [];
                    cleanMEG.trialinfo.triggerchannel([3,10,11,12,13,20,22,25,26],:) = [];
                    cleanMEG.trialinfo.responsechannel([3,10,11,12,13,20,22,25,26],:) = [];
                    cleanMEG.trialinfo.triggerlabel([3,10,11,12,13,20,22,25,26]) = [];
                    cleanMEG.trialinfo.response([3,10,11,12,13,20,22,25,26]) = [];
                    cleanMEG.trialinfo.response_label([3,10,11,12,13,20,22,25,26]) = [];
                    cleanMEG.trialinfo.balldesign([3,10,11,12,13,20,22,25,26]) = [];
                    cleanMEG.trialinfo.balldesign_short([3,10,11,12,13,20,22,25,26]) = [];
                 elseif 1 == strcmp(files(ii).name, 'Neu_Dislike500_3.mat')
                    cleanMEG.trial([11]) = [];
                    cleanMEG.time([11]) = [];
                    cleanMEG.trialinfo.triggerchannel([11],:) = [];
                    cleanMEG.trialinfo.responsechannel([11],:) = [];
                    cleanMEG.trialinfo.triggerlabel([11]) = [];
                    cleanMEG.trialinfo.response([11]) = [];
                    cleanMEG.trialinfo.response_label([11]) = [];
                    cleanMEG.trialinfo.balldesign([11]) = [];
                    cleanMEG.trialinfo.balldesign_short([11]) = [];    
                    
                 elseif 1 == strcmp(files(ii).name, 'Neu_Like500_1.mat')
                     % 39 trials, davon keine trials entfernt
                     
                     
                     
                 elseif 1 == strcmp(files(ii).name, 'Neu_Like500_2.mat')    
                   % unklar, welche ich hier entfernt habe
                    
                elseif 1 == strcmp(files(ii).name, 'Neu_Like500_3.mat')
                    cleanMEG.trial([15]) = [];
                    cleanMEG.time(15) = [];
                    cleanMEG.trialinfo.triggerchannel(15,:) = [];
                    cleanMEG.trialinfo.responsechannel(15,:) = [];
                    cleanMEG.trialinfo.triggerlabel(15) = [];
                    cleanMEG.trialinfo.response(15) = [] ;
                    cleanMEG.trialinfo.response_label(15) = [];
                    cleanMEG.trialinfo.balldesign(15) = [];
                    cleanMEG.trialinfo.balldesign_short(15) = [];
                    
                end
                
                cleanMEG.additional_cleaning = 'yes';
            end



             %% plot and save avg after cleaning
            cfgn                = [];
            cfgn.parameter      = 'trial';
            cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
            cfgn.vartrllength   = 2;

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

