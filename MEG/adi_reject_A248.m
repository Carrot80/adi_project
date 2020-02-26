

function adi_rejectvisual_MEG_extra(path2cleanfile, subject)


PathFigure_after_cleaning = [path2cleanfile, 'figures\after_cleaning\'];   

files   =   dir(fullfile(path2cleanfile, '*.mat'));
size_files = size(files);

for ii = 1:(size_files(1,1))

load ([files(ii).folder, filesep, files(ii).name], 'cleanMEG')
             
bad_chan = {'A248'};
bad_chan_ind = [];
for kk = 1:length(bad_chan)
     bad_chan_ind(kk) = find(strcmp(cleanMEG.label, bad_chan{kk}));
end


if ~isnan(cleanMEG.trial{kk}(bad_chan_ind,:))
    for kk = 1:length(cleanMEG.trial)
        cleanMEG.trial{kk}(bad_chan_ind,:) = NaN;
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

            saveas (gcf, [PathFigure_after_cleaning files(ii).name(1:end-4) '_A248_rejected'], 'fig') % save EEG after cleaning figure to directory     
            save ([path2cleanfile files(ii).name], 'cleanMEG'); % save mat-file after cleaning    
            clear cleanMEG
            disp(['A248 additionally rejected in subject ' subject ])
    
else
    return
end

end

        
               
                
            
        
          
       

    


end

