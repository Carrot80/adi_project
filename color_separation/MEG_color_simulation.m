%% trials sollen gruppiert werden anhand der Muster bzw. Farben der Bälle
% adi add trigger codes 2 cleanData, nachträglich: 

clear
path2data = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2data);
subject_list([1 2],:) = [];
subject_list([1 2 5 9 14],:) = [];

color.red_white = {'rwv', 'rwf', 'rws'};
color.grey_green = {'ggv', 'ggf', 'ggs'};
color.yellow_blue = {'gbv', 'gbf', 'gbs'};

triggercode_labels = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
triggercodes = {'102', '104', '106', '108', '4196', '4198', '4200', '4202', '4204'};

add_triggercode_labels2data(path2data, subject_list, color, triggercodes, triggercode_labels, num2str(1));
add_triggercode_labels2data(path2data, subject_list, color, triggercodes, triggercode_labels, num2str(2));
add_triggercode_labels2data(path2data, subject_list, color, triggercodes, triggercode_labels, num2str(3));
%%
% adi add trigger codes 2 interpolated data, nachträglich: 

clear
path2data = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
path_add_triggers = 'MEG_analysis\noisereduced\1_95Hz\02_interpolated';
path_trigger_sources = 'MEG_analysis\noisereduced\1_95Hz\01_clean\';
subject_list = dir (path2data);
subject_list([1 2],:) = [];
% subject_list([1 2 5 9 14],:) = [];

color.red_white = {'rwv', 'rwf', 'rws'};
color.grey_green = {'ggv', 'ggf', 'ggs'};
color.yellow_blue = {'gbv', 'gbf', 'gbs'};

triggercode_labels = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
triggercodes = {'102', '104', '106', '108', '4196', '4198', '4200', '4202', '4204'};

add_triggercode_labels2interpolated_data(path2data, subject_list, path_add_triggers, 'cleanMEG_interp', path_trigger_sources, 'cleanMEG', color, triggercodes, triggercode_labels);

%%
    clear   
%     brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_95Hz';
    
    trigger.pattern        = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
    trigger.triggerchannel = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
    trigger.eprime        = [102, 104, 106, 108, 101, 103, 105, 107, 109];

    
    for i = 1: length(ListSubj)  
        path2cleanfile = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\' filter '\01_clean\']);
        pathInterpolated = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\' filter '\02_interpolated\']);
        
        adi_add_trialinfo2subjects (path2cleanfile, pathInterpolated, trigger)
    end
 

%% group trials of each run according 2 three patterns volley, space, soccerball 

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
freq = 'bp1-45Hz';
load('U:\My Documents\MATLAB\atlas_source_indices.mat');
condition.Volley = 1;
condition.Space = 2;
condition.Soccer = 3;

virtsens_all_subj = struct('trial', [], 'time', [], 'balldesign', [], 'ballpattern', [], 'trialnumbers', [], 'triggers', [], 'label', [], 'fsample', []);

for i=1:7%length(subject_list)
    
    vs_allRuns = struct();
    path2sensordata = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
    path2spatialfilter = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\spatialfilter\'];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\pattern_comparison'];
    
%     if ~exist   ([dir_out '\virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'file')
    [vs_allRuns] = adi_group_patterns(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(1));
    [vs_allRuns] = adi_group_patterns(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(2));
    [vs_allRuns] = adi_group_patterns(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(3));

    run = fields(vs_allRuns);

    switch size(fields(vs_allRuns),1)
        case 3 
            vs_allRuns_appended.trial = [vs_allRuns.run1.trial vs_allRuns.run2.trial vs_allRuns.run3.trial];
            vs_allRuns_appended.time =  [vs_allRuns.run1.time vs_allRuns.run2.time vs_allRuns.run3.time];
            vs_allRuns_appended.balldesign =  [vs_allRuns.run1.balldesign vs_allRuns.run2.balldesign vs_allRuns.run3.balldesign];
            vs_allRuns_appended.ballpattern =  [vs_allRuns.run1.ballpattern vs_allRuns.run2.ballpattern vs_allRuns.run3.ballpattern];
            vs_allRuns_appended.trialnumbers =  [vs_allRuns.run1.trialnumbers vs_allRuns.run2.trialnumbers vs_allRuns.run3.trialnumbers];
            vs_allRuns_appended.triggers =  [vs_allRuns.run1.triggers vs_allRuns.run2.triggers vs_allRuns.run3.triggers];         
        case 2
            vs_allRuns_appended.trial = [vs_allRuns.run1.trial vs_allRuns.run2.trial];
            vs_allRuns_appended.time =  [vs_allRuns.run1.time vs_allRuns.run2.time];
            vs_allRuns_appended.balldesign =  [vs_allRuns.run1.balldesign vs_allRuns.run2.balldesign];
            vs_allRuns_appended.ballpattern =  [vs_allRuns.run1.ballpattern vs_allRuns.run2.ballpattern];
            vs_allRuns_appended.trialnumbers =  [vs_allRuns.run1.trialnumbers vs_allRuns.run2.trialnumbers];
            vs_allRuns_appended.triggers =  [vs_allRuns.run1.triggers vs_allRuns.run2.triggers];         

    end
    
    vs_allRuns_appended.label =  vs_allRuns.run1.label;
    vs_allRuns_appended.fsample =  vs_allRuns.run1.fsample;
    clearvars vs_allRuns
    % sanity check Nr. 3:
   
    cfg = [];
    avg = ft_timelockanalysis(cfg, vs_allRuns_appended);
    figure;
    plot(avg.time, avg.avg)
    savefig([dir_out filesep 'sanity_check\' 'avg_virtsens_allRuns_' freq '.fig'])
    close
    
    pathAppended = [dir_out '\virtsens_ns\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end
   
    save ([pathAppended 'virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended', '-v7.3');

%     else
%         load([dir_out '\virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended')
%     end
    
    virtsens_all_subj(i) = vs_allRuns_appended;
    clearvars vs_allRuns_appended
end

    
   session = virtsens_all_subj(1);
   for k =2:length(virtsens_all_subj)
       session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
       session.time = cat(2,session.time, virtsens_all_subj(k).time);
       session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
       session.ballpattern = cat(2,session.ballpattern, virtsens_all_subj(k).ballpattern);
       session.triggers = cat(2,session.triggers, virtsens_all_subj(k).triggers);
       session.trialnumbers = cat(2,session.trialnumbers, virtsens_all_subj(k).trialnumbers);
   end
    
   clearvars virtsens_all_subj
   
for i = 1:length(session.trial)
    temp = session.trial{1,i};
    session.data(i,:,:) = temp;
    clear temp
end
   session = rmfield(session, 'trial');
   session.label = session.ballpattern;
    
    group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\pattern_comparison\guggenmos_svm_results\';
    svm_guggenmos_pattern (session, group_path, freq, condition, atlas_downsampled)
    
    
    %%
    
    
    
    
    
    
 