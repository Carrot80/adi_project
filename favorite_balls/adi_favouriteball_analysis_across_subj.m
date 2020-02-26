function [] = adi_favouriteball_analysis(subject_dir, path2file)

% list_favorite_balls = import_favouriteballlist('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\list_favourite_balls.txt');
data_favorite1 =  struct('trial', [], 'time', [], 'balldesign_short', [], 'labels', [],  'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_favorite2 =  struct('trial', [], 'time', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_favorite3 =  struct('trial', [], 'time', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_favorite4 =  struct('trial', [], 'time', [], 'balldesign_short', [], 'labels', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);
data_dislike = struct('trial', [], 'time', [], 'balldesign_short', [], 'label', [], 'fsample', [] ,'grad', [], 'cfg', [], 'subject', []);

for ii = 1:length(subject_dir)

    load ([subject_dir(ii).folder filesep subject_dir(ii).name '\favorite_balldesign_1.mat'])   
    load ([subject_dir(ii).folder filesep subject_dir(ii).name '\favorite_balldesign_2.mat'])   
    load ([subject_dir(ii).folder filesep subject_dir(ii).name '\favorite_balldesign_3.mat'])   
    load ([subject_dir(ii).folder filesep subject_dir(ii).name '\favorite_balldesign_4.mat'])   
    load ([subject_dir(ii).folder filesep subject_dir(ii).name '\session_dislike.mat'])   
    
    data_favorite1(ii) = session_favorite1;
    data_favorite2(ii) = session_favorite2;
    data_favorite3(ii) = session_favorite3;
    data_favorite4(ii) = session_favorite4;
    data_dislike(ii) = session_dislike;
    
%     figure;
%     plot(session_dislike.time{1,1}, rms_favorite1_data)
%     hold on; plot(session.time{1,1}, rms_favorite2_data)
%     hold on; plot(session.time{1,1}, rms_favorite3_data)
%     hold on; plot(session.time{1,1}, rms_favorite4_data)
%     hold on; plot(session.time{1,1}, rms_dislike_data)
%     legend({'fav1', 'fav2', 'fav3', 'fav4', 'unpreferred'})
%     
%     fav_1to4 = [favorite1_data; favorite2_data; favorite3_data; favorite4_data];
%     
%     figure;
%     plot(session.time{1,1}, rms(squeeze(mean(fav_1to4))))
%     hold on; plot(session.time{1,1}, rms_dislike_data)
    
%     %% separate unpreferred balls:
%     
%     balldesigns_unpreferred = unique(session_dislike.balldesign_short);
%     
%     for kk = 1:length(balldesigns_unpreferred)
%         dislike.(balldesigns_unpreferred{kk}).ind = find(strcmp(session_dislike.balldesign_short, balldesigns_unpreferred{kk}));
%         dislike.(balldesigns_unpreferred{kk}).data = kh_trial2dat(session_dislike.trial(dislike.(balldesigns_unpreferred{kk}).ind));
%     end
%        
%     figure; hold on;
%     for kk = 1:length(balldesigns_unpreferred)
%         plot(session.time{1,1}, rms(squeeze(mean(dislike.(balldesigns_unpreferred{kk}).data))))
%     end
%    
%     hold on; plot(session.time{1,1}, rms_favorite1_data)
%     hold on; plot(session.time{1,1}, rms_favorite2_data)
%     hold on; plot(session.time{1,1}, rms_favorite3_data)
%     hold on; plot(session.time{1,1}, rms_favorite4_data)
%     
    
%     hold on; plot(session.time{1,1}, rms_dislike_data)
%     
%     plot(session.time{1,1}, rms(squeeze(mean(fav_1to4))))
%     
    
%     figure;
%     plot(session.time{1,1}, rms(squeeze(mean(fav_1to4))))
%     hold on; plot(session.time{1,1}, rms_dislike_data)
    
%     close all
      
     
end

dat_favorite1 = data_favorite1(1);

for kk = 2:length(data_favorite1)
    dat_favorite1.trial = cat(2, dat_favorite1.trial, data_favorite1(kk).trial);
    dat_favorite1.time = cat(2, dat_favorite1.time, data_favorite1(kk).time);
    dat_favorite1.balldesign_short = cat(2, dat_favorite1.balldesign_short, data_favorite1(kk).balldesign_short);
    dat_favorite1.subject = cat(2, dat_favorite1.subject, data_favorite1(kk).subject);
    dat_favorite1.labels = cat(2, dat_favorite1.labels, data_favorite1(kk).labels);
end

dat_favorite2 = data_favorite2(1);

for kk = 2:length(data_favorite2)
    dat_favorite2.trial = cat(2, dat_favorite2.trial, data_favorite2(kk).trial);
    dat_favorite2.time = cat(2, dat_favorite2.time, data_favorite2(kk).time);
    dat_favorite2.balldesign_short = cat(2, dat_favorite2.balldesign_short, data_favorite2(kk).balldesign_short);
    dat_favorite2.subject = cat(2, dat_favorite2.subject, data_favorite2(kk).subject);
    dat_favorite2.labels = cat(2, dat_favorite2.labels, data_favorite2(kk).labels);
end

dat_favorite3 = data_favorite3(1);

for kk = 2:length(data_favorite3)
    dat_favorite3.trial = cat(2, dat_favorite3.trial, data_favorite3(kk).trial);
    dat_favorite3.time = cat(2, dat_favorite3.time, data_favorite3(kk).time);
    dat_favorite3.balldesign_short = cat(2, dat_favorite3.balldesign_short, data_favorite3(kk).balldesign_short);
    dat_favorite3.subject = cat(2, dat_favorite3.subject, data_favorite3(kk).subject);
    dat_favorite3.labels = cat(2, dat_favorite3.labels, data_favorite3(kk).labels);
end

dat_favorite4 = data_favorite4(1);

for kk = 2:length(data_favorite4)
    dat_favorite4.trial = cat(2, dat_favorite4.trial, data_favorite4(kk).trial);
    dat_favorite4.time = cat(2, dat_favorite4.time, data_favorite4(kk).time);
    dat_favorite4.balldesign_short = cat(2, dat_favorite4.balldesign_short, data_favorite4(kk).balldesign_short);
    dat_favorite4.subject = cat(2, dat_favorite4.subject, data_favorite4(kk).subject);
    dat_favorite4.labels = cat(2, dat_favorite4.labels, data_favorite4(kk).labels);
end

dat_dislike = data_dislike(1);

for kk = 2:length(data_dislike)
    dat_dislike.trial = cat(2, dat_dislike.trial, data_dislike(kk).trial);
    dat_dislike.time = cat(2, dat_dislike.time, data_dislike(kk).time);
    dat_dislike.balldesign_short = cat(2, dat_dislike.balldesign_short, data_dislike(kk).balldesign_short);
    dat_dislike.subject = cat(2, dat_dislike.subject, data_dislike(kk).subject);
%     dat_dislike.labels = cat(2, dat_dislike.labels, data_dislike(kk).labels);
end

rand_ = randperm(numel(dat_dislike.trial), numel(dat_favorite1.trial));
 

dat_dislike.trial = dat_dislike.trial(rand_); 
dat_dislike.time = dat_dislike.time(rand_);

cfg = [];
avg_favorit_1 = ft_timelockanalysis(cfg, dat_favorite1);
avg_favorit_2 = ft_timelockanalysis(cfg, dat_favorite2);
avg_favorit_3 = ft_timelockanalysis(cfg, dat_favorite3);
avg_favorit_4 = ft_timelockanalysis(cfg, dat_favorite4);
avg_dislike = ft_timelockanalysis(cfg, dat_dislike);



figure
plot(avg_favorit_1.time, rms(avg_favorit_1.avg))
hold on
plot(avg_favorit_1.time, rms(avg_favorit_2.avg))
hold on
plot(avg_favorit_1.time, rms(avg_favorit_3.avg))
hold on
plot(avg_favorit_1.time, rms(avg_favorit_4.avg))
hold on
plot(avg_favorit_1.time, rms(avg_dislike.avg))
box off
legend({'favorite rank 1'; 'favorite rank 2'; 'favorite rank 3'; 'favorite rank 4'; 'unfavored'})
set(legend,'Box','off')


figure; 
plot(avg_favorit_1.time, rms(avg_favorit_1.avg))
hold on; plot(avg_favorit_1.time, rms(avg_dislike.avg))
box off
legend({'favorite rank 1'; 'unfavored'})
set(legend,'Box','off')

figure
plot(avg_favorit_1.time, rms(avg_favorit_2.avg))
hold on; plot(avg_favorit_1.time, rms(avg_dislike.avg))
box off
legend({'favorite rank 2'; 'unfavored'})
set(legend,'Box','off')

figure
plot(avg_favorit_1.time, rms(avg_favorit_3.avg))
hold on; plot(avg_favorit_1.time, rms(avg_dislike.avg))
box off
legend({'favorite rank 3'; 'unfavored'})
set(legend,'Box','off')

figure
plot(avg_favorit_1.time, rms(avg_favorit_4.avg))
hold on; plot(avg_favorit_1.time, rms(avg_dislike.avg))
box off
legend({'favorite rank 4'; 'unfavored'})
set(legend,'Box','off')

%%

dat_favorite_all = dat_favorite1;
dat_favorite_all.trial = [dat_favorite1.trial dat_favorite2.trial dat_favorite3.trial dat_favorite4.trial];
dat_favorite_all.time = [dat_favorite1.time dat_favorite2.time dat_favorite3.time dat_favorite4.time];

dat_dislike_all = data_dislike(1);

for kk = 2:length(data_dislike)
    dat_dislike_all.trial = cat(2, dat_dislike_all.trial, data_dislike(kk).trial);
    dat_dislike_all.time = cat(2, dat_dislike_all.time, data_dislike(kk).time);
    dat_dislike_all.balldesign_short = cat(2, dat_dislike_all.balldesign_short, data_dislike(kk).balldesign_short);
    dat_dislike_all.subject = cat(2, dat_dislike_all.subject, data_dislike(kk).subject);
%     dat_dislike.labels = cat(2, dat_dislike.labels, data_dislike(kk).labels);
end

rand_ = randperm(numel(dat_dislike_all.trial), numel(dat_favorite_all.trial));

dat_dislike_all.trial = dat_dislike_all.trial(rand_); 
dat_dislike_all.time = dat_dislike_all.time(rand_);

cfg = [];
avg_favorit_all = ft_timelockanalysis(cfg, dat_favorite_all);
avg_dislike_all = ft_timelockanalysis(cfg, dat_dislike_all);

figure;
plot(avg_favorit_1.time, rms(avg_favorit_all.avg))
hold on
plot(avg_favorit_1.time, rms(avg_dislike_all.avg))
legend({'favorites_all'; 'dislike'})



%% favorite 2-4 zusammenfassen:
avg_favorit2_4.avg = mean(cat(3, avg_favorit_2.avg, avg_favorit_3.avg, avg_favorit_4.avg),3);
figure; hold on; box off
plot(avg_favorit_1.time, rms(avg_favorit_1.avg))
plot(avg_favorit_1.time, rms(avg_favorit2_4.avg))
plot(avg_favorit_1.time, rms(avg_dislike.avg), 'g')
legend({'favorite rank 1';  'favorite rank 2 - 4'; 'unfavored'})
set(legend,'Box','off')
xlabel('time')
ylabel('')

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

