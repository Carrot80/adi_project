function [] = adi_favouriteball_analysis(subject_dir, path2file)

list_favorite_balls = import_favouriteballlist('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\list_favourite_balls.txt');
data_favorite1 =  struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [], 'labels', [],  'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_favorite2 =  struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_favorite3 =  struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_favorite4 =  struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_dislike = struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);

for ii = 1:2%:length(subject_dir)

    dir_data = dir([subject_dir(ii).folder filesep subject_dir(ii).name filesep path2file '*.mat']);

    for kk = 1:length(dir_data)
        load ([dir_data(kk).folder filesep dir_data(kk).name], 'cleanMEG_interp')
        for pp = 1:length(cleanMEG_interp.trial)
            cleanMEG_interp.trialinfo.run{pp} = dir_data(kk).name(end-4);
        end
        if isfield(cleanMEG_interp, 'additional_cleaning')
            cleanMEG_interp = rmfield(cleanMEG_interp, 'additional_cleaning');
        end
        data(kk) = cleanMEG_interp;
        clear cleanMEG_interp   
    end
    
    session = data(1);
    session.response_label = data(1).trialinfo.response_label;
    session.run = data(1).trialinfo.run;
    session.balldesign = data(1).trialinfo.balldesign_short;
    session = rmfield(session, 'trialinfo');

    for kk = 2:length(data)
        session.trial = cat(2, session.trial, data(kk).trial);
        session.time = cat(2, session.time, data(kk).time);
        session.response_label = cat(2, session.response_label, data(kk).trialinfo.response_label);
        session.balldesign = cat(2, session.balldesign, data(kk).trialinfo.balldesign_short);
        session.run = cat(2, session.run, data(kk).trialinfo.run);
    end

    clear data
    for kk = 1:length(session.response_label)
        switch session.response_label{kk}
            case 'like' % 'Volley'
                session.labels(kk) = 1;
            case 'Neu_Like' % 'Volley'
                session.labels(kk) = 1;
            case 'dislike' % 'Space'
                session.labels(kk) = 2;
            case 'Neu_Dislike' % 'Space'
                session.labels(kk) = 2;
        end
    end    

   [session] = adi_bpfilter(session, 'bp1_45Hz');     
   
    data_trials = kh_trial2dat(session.trial);
    [FourRef,Fref]=fftBasic(squeeze(mean(data_trials,1)), session.fsample);
    figure;
    plot(FourRef, Fref)
   
  [session] = adi_ztrans_sensorspace(session);   
   
   [session] = recode_favorite_balldesigns(session, list_favorite_balls, subject_dir(ii).name);
   
    indx_ = [];
    for pp = 1:length(list_favorite_balls)
        indx_(pp) = strcmp(list_favorite_balls(pp,1), subject_dir(ii).name);
    end

    indx_subject = find(indx_);

    favorite_ball_1 = list_favorite_balls(indx_subject,2);
    favorite_ball_2 = list_favorite_balls(indx_subject,3);
    favorite_ball_3 = list_favorite_balls(indx_subject,4);
    favorite_ball_4 = list_favorite_balls(indx_subject,5);
    
    trls_favorite_1 = find(strcmp(session.balldesign, favorite_ball_1));
    trls_favorite_2 = find(strcmp(session.balldesign, favorite_ball_2));
    trls_favorite_3 = find(strcmp(session.balldesign, favorite_ball_3));
    trls_favorite_4 = find(strcmp(session.balldesign, favorite_ball_4));
    
    
    
    
    %% session_favorite1
    session_favorite1.trial = session.trial(trls_favorite_1);
    session_favorite1.time = session.time(trls_favorite_1);
    session_favorite1.balldesign_short = session.balldesign(trls_favorite_1);
    session_favorite1.labels = session.labels(trls_favorite_1);
    session_favorite1.label = session.label;
    
    session_favorite1.fsample = session.fsample;
    session_favorite1.grad = session.grad;
    session_favorite1.cfg = session.cfg;
    
    session_favorite1 = setfield(session_favorite1, 'subject', []);
    for pp = 1:length(trls_favorite_1)
        session_favorite1.subject{pp} = subject_dir(ii).name;
    end
    clear favorite_ball_1
    
      %% session_favorite2
    session_favorite2.trial = session.trial(trls_favorite_2);
    session_favorite2.time = session.time(trls_favorite_2);
    session_favorite2.balldesign_short = session.balldesign(trls_favorite_2);
    session_favorite2.labels = session.labels(trls_favorite_2);
    session_favorite2.label = session.label;
    
    session_favorite2.fsample = session.fsample;
    session_favorite2.grad = session.grad;
    session_favorite2.cfg = session.cfg;
    
    session_favorite2 = setfield(session_favorite2, 'subject', []);
    for pp = 1:length(trls_favorite_2)
        session_favorite2.subject{pp} = subject_dir(ii).name;
    end
    clear favorite_ball_2
    
       %% session_favorite3
    session_favorite3.trial = session.trial(trls_favorite_3);
    session_favorite3.time = session.time(trls_favorite_3);
    session_favorite3.balldesign_short = session.balldesign(trls_favorite_3);
    session_favorite3.labels = session.labels(trls_favorite_3);
    session_favorite3.label = session.label;
    session_favorite3.fsample = session.fsample;
    session_favorite3.grad = session.grad;
    session_favorite3.cfg = session.cfg;
    
    session_favorite3 = setfield(session_favorite3, 'subject', []);
    for pp = 1:length(trls_favorite_3)
        session_favorite3.subject{pp} = subject_dir(ii).name;
    end
    clear favorite_ball_3
    
      %% session_favorite4
    session_favorite4.trial = session.trial(trls_favorite_4);
    session_favorite4.time = session.time(trls_favorite_4);
    session_favorite4.balldesign_short = session.balldesign(trls_favorite_4);
    session_favorite4.labels = session.labels(trls_favorite_4);
    session_favorite4.label = session.label;
    session_favorite4.fsample = session.fsample;
    session_favorite4.grad = session.grad;
    session_favorite4.cfg = session.cfg;
    
    session_favorite4 = setfield(session_favorite4, 'subject', []);
    for pp = 1:length(trls_favorite_4)
        session_favorite4.subject{pp} = subject_dir(ii).name;
    end
    
    clear favorite_ball_4
    
    %% session_dislike:
    
    session_dislike = [];
    session_dislike.trial = session.trial;
    session_dislike.trial([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]) = [];
    session_dislike.time = session.time;
    session_dislike.time([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]) = [];
    session_dislike.balldesign_short = session.balldesign;
    session_dislike.balldesign_short([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]) = [];
    session_dislike.label = session.label;
    session_dislike.fsample = session.fsample;
    session_dislike.grad = session.grad;
    session_dislike.cfg = session.cfg;
    
    session_dislike = setfield(session_dislike, 'subject', []);
    for pp = 1:length(session_dislike.trial)
        session_dislike.subject{pp} = subject_dir(ii).name;
    end
    
    favorite1_data = kh_trial2dat(session_favorite1.trial);
    favorite2_data = kh_trial2dat(session_favorite2.trial);
    favorite3_data = kh_trial2dat(session_favorite3.trial);
    favorite4_data = kh_trial2dat(session_favorite4.trial);
    dislike_data = kh_trial2dat(session_dislike.trial);
    
    rms_favorite1_data = rms(squeeze(mean(favorite1_data,1)));
    rms_favorite2_data = rms(squeeze(mean(favorite2_data,1))); 
    rms_favorite3_data = rms(squeeze(mean(favorite3_data,1))); 
    rms_favorite4_data = rms(squeeze(mean(favorite4_data,1))); 
    rms_dislike_data = rms(squeeze(mean(dislike_data,1))); 
    
    figure;
    plot(session.time{1,1}, rms_favorite1_data)
    hold on; plot(session.time{1,1}, rms_favorite2_data)
    hold on; plot(session.time{1,1}, rms_favorite3_data)
    hold on; plot(session.time{1,1}, rms_favorite4_data)
    hold on; plot(session.time{1,1}, rms_dislike_data)
    legend({'fav1', 'fav2', 'fav3', 'fav4', 'unpreferred'})
    
    fav_1to4 = [favorite1_data; favorite2_data; favorite3_data; favorite4_data];
    
    figure;
    plot(session.time{1,1}, rms(squeeze(mean(fav_1to4))))
    hold on; plot(session.time{1,1}, rms_dislike_data)
    
    %% separate unpreferred balls:
    
    balldesigns_unpreferred = unique(session_dislike.balldesign_short);
    
    for kk = 1:length(balldesigns_unpreferred)
        dislike.(balldesigns_unpreferred{kk}).ind = find(strcmp(session_dislike.balldesign_short, balldesigns_unpreferred{kk}));
        dislike.(balldesigns_unpreferred{kk}).data = kh_trial2dat(session_dislike.trial(dislike.(balldesigns_unpreferred{kk}).ind));
    end
       
    figure; hold on;
    for kk = 1:length(balldesigns_unpreferred)
        plot(session.time{1,1}, rms(squeeze(mean(dislike.(balldesigns_unpreferred{kk}).data))))
    end
   
    hold on; plot(session.time{1,1}, rms_favorite1_data)
    hold on; plot(session.time{1,1}, rms_favorite2_data)
    hold on; plot(session.time{1,1}, rms_favorite3_data)
    hold on; plot(session.time{1,1}, rms_favorite4_data)
    
    
%     hold on; plot(session.time{1,1}, rms_dislike_data)
%     
%     plot(session.time{1,1}, rms(squeeze(mean(fav_1to4))))
%     
    
%     figure;
%     plot(session.time{1,1}, rms(squeeze(mean(fav_1to4))))
%     hold on; plot(session.time{1,1}, rms_dislike_data)
    
%     data_favorite1(ii) = session_favorite1;
%     data_favorite2(ii) = session_favorite2;
%     data_favorite3(ii) = session_favorite3;
%     data_favorite4(ii) = session_favorite4;
%     data_dislike(ii) = session_dislike;
    
    
%     close all
      
    if ~exist(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name filesep])
        mkdir((['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name])) 
    end
    save(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name filesep 'favorite_balldesign_1.mat'], 'session_favorite1')
    save(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name filesep 'favorite_balldesign_2.mat'], 'session_favorite2')
    save(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name filesep 'favorite_balldesign_3.mat'], 'session_favorite3')
    save(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name filesep 'favorite_balldesign_4.mat'], 'session_favorite4')
    save(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\' subject_dir(ii).name filesep 'session_dislike.mat'], 'session_dislike')
 
    clear dir_data data session session_favorite1 session_favorite2 session_favorite3 session_favorite4 session_dislike
    clear trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4 
    
end

dat_favorite1 = data_favorite1(1);

for kk = 2:length(data_favorite1)
    dat_favorite1.trial = cat(2, dat_favorite1.trial, data_favorite1(kk).trial);
    dat_favorite1.time = cat(2, dat_favorite1.time, data_favorite1(kk).time);
    dat_favorite1.response_label = cat(2, dat_favorite1.response_label, data_favorite1(kk).response_label);
    dat_favorite1.balldesign_short = cat(2, dat_favorite1.balldesign_short, data_favorite1(kk).balldesign_short);
    dat_favorite1.subject = cat(2, dat_favorite1.subject, data_favorite1(kk).subject);
     dat_favorite1.labels = cat(2, dat_favorite1.labels, data_favorite1(kk).labels);
end

dat_favorite2 = data_favorite2(1);

for kk = 2:length(data_favorite2)
    dat_favorite2.trial = cat(2, dat_favorite2.trial, data_favorite2(kk).trial);
    dat_favorite2.time = cat(2, dat_favorite2.time, data_favorite2(kk).time);
    dat_favorite2.response_label = cat(2, dat_favorite2.response_label, data_favorite2(kk).response_label);
    dat_favorite2.balldesign_short = cat(2, dat_favorite2.balldesign_short, data_favorite2(kk).balldesign_short);
    dat_favorite2.subject = cat(2, dat_favorite2.subject, data_favorite2(kk).subject);
    dat_favorite2.labels = cat(2, dat_favorite2.labels, data_favorite2(kk).labels);
end

dat_favorite3 = data_favorite3(1);

for kk = 2:length(data_favorite3)
    dat_favorite3.trial = cat(2, dat_favorite3.trial, data_favorite3(kk).trial);
    dat_favorite3.time = cat(2, dat_favorite3.time, data_favorite3(kk).time);
    dat_favorite3.response_label = cat(2, dat_favorite3.response_label, data_favorite3(kk).response_label);
    dat_favorite3.balldesign_short = cat(2, dat_favorite3.balldesign_short, data_favorite3(kk).balldesign_short);
    dat_favorite3.subject = cat(2, dat_favorite3.subject, data_favorite3.labels, data_favorite3(kk).labels);
end

dat_favorite4 = data_favorite4(1);

for kk = 2:length(data_favorite4)
    dat_favorite4.trial = cat(2, dat_favorite4.trial, data_favorite4(kk).trial);
    dat_favorite4.time = cat(2, dat_favorite4.time, data_favorite4(kk).time);
    dat_favorite4.response_label = cat(2, dat_favorite3.response_label, data_favorite4(kk).response_label);
    dat_favorite4.balldesign_short = cat(2, dat_favorite4.balldesign_short, data_favorite4(kk).balldesign_short);
    dat_favorite4.subject = cat(2, dat_favorite4.subject, data_favorite4.labels, data_favorite4(kk).labels);
end

cfg=[];
avg_favorit_1 = ft_timelockanalysis(cfg, dat_favorite1)
avg_favorit_2 = ft_timelockanalysis(cfg, dat_favorite2)
avg_favorit_3 = ft_timelockanalysis(cfg, dat_favorite3)
avg_favorit_4 = ft_timelockanalysis(cfg, dat_favorite4)


figure
plot(avg_favorit_1.time, rms(avg_favorit_1.avg))
hold on
plot(avg_favorit_1.time, rms(avg_favorit_2.avg))
hold on
plot(avg_favorit_1.time, rms(avg_favorit_3.avg))
hold on
plot(avg_favorit_1.time, rms(avg_favorit_4.avg))
legend({'favorite 1'; 'favorite 2'; 'favorite 3'; 'favorite 4'})

end

function [session] = recode_favorite_balldesigns(session, list_favorite_balls, subject)

    indx_ = [];
    for pp = 1:length(list_favorite_balls)
        indx_(pp) = strcmp(list_favorite_balls(pp,1), subject);
    end

    indx_subject = find(indx_);
    
    favorite_ball_1 = list_favorite_balls(indx_subject,2);
    favorite_ball_2 = list_favorite_balls(indx_subject,3);
    favorite_ball_3 = list_favorite_balls(indx_subject,4);
    favorite_ball_4 = list_favorite_balls(indx_subject,5);
    
    for pp = 1:numel(session.balldesign)
        session.balldesign(pp) = session.balldesign{pp};
    end
    trls_favorite_1 = find(strcmp(session.balldesign, favorite_ball_1));
    trls_favorite_2 = find(strcmp(session.balldesign, favorite_ball_2));
    trls_favorite_3 = find(strcmp(session.balldesign, favorite_ball_3));
    trls_favorite_4 = find(strcmp(session.balldesign, favorite_ball_4));
    
    trials_like = session.trial([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]);
    trials_dislike = session.trial;
    trials_dislike([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]) = [];
    balldesign_dislike = session.balldesign;
    balldesign_dislike([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]) = [];
    balldesign_like = session.balldesign([trls_favorite_1 trls_favorite_2 trls_favorite_3 trls_favorite_4]);
    
    session.trial = [];
    session.trial = [trials_like trials_dislike];
    session.labels = [ones(1, numel(trials_like)) 2*ones(1, numel(trials_dislike))];
    session.response_label = [];
    session.balldesign = [];
    session.balldesign = [balldesign_like balldesign_dislike];
    session.run = [];
    
end


function [data_bpfreq_res_sel] = adi_bpfilter(filename, bpname)


switch bpname
    case 'bp1_95Hz'
       data_bpfreq_res_sel =  filename;
        return
    case 'bp1_45Hz'
        bpfreq = [1 45];
    case'1_5_45Hz'
        bpfreq = [1.5 45];
    case 'bp2_45Hz'
        bpfreq = [2 45];
    case 'bp3_45Hz'
        bpfreq = [3 45];
    case 'delta'
        bpfreq = 4;
    case 'theta'
        bpfreq = [4 8];
    case 'alpha'
        bpfreq = [8 13];
    case 'beta'
        bpfreq = [13 25];
    case 'low_gamma'
        bpfreq = [25 45];
    case 'high_gamma'
        bpfreq = [55 90];
    case 'bp10-45Hz'
        bpfreq = [10 45];
end

cfg = [];
cfg.keeptrials = 'yes';
cfg.vartrllength = 2;

cfg.trials  = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') %|| 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
    cfg.demean = 'yes';
    cfg.baselinewindow  = [-0.5 -0.030];
end

try
    [data_bpfreq] = ft_preprocessing(cfg, filename); 
    [warnMsg, warnID] = lastwarn;
    if ~isempty(warnMsg)
       warnMsg
    end
catch
    
    warnMsg
end

cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);

data_bpfreq_res_sel.grad.label(249:end) = [];
data_bpfreq_res_sel.grad.chanori(249:end, :) = [];
data_bpfreq_res_sel.grad.chanpos(249:end, :) = [];
data_bpfreq_res_sel.grad.tra(249:end, :) = [];
data_bpfreq_res_sel.label(249:end) = [];

fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_res_sel);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);

for k=1:length(diff_fieldnames)
    data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end

if isfield(data_bpfreq_res_sel, 'sampleinfo')
    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'sampleinfo');
end

% fn_filename{end+1}='cfg';
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, fn_filename);

if ~isfield(data_bpfreq_res_sel, 'additional_cleaning')
    data_bpfreq_res_sel = setfield(data_bpfreq_res_sel, 'additional_cleaning', 'no');
end

right_fieldorder = {'trial'; 'time'; 'label'; 'balldesign'; 'labels'; 'response_label'; 'ChannelFlag_Bst'; 'additional_cleaning'; 'cfg'; 'dimord'; 'fsample'; 'grad'; 'run'};
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, right_fieldorder);

clearvars filename data_bpfreq data_bpfreq_res 

 
end

